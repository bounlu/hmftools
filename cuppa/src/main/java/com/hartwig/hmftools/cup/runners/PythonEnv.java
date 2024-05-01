package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;

public class PythonEnv
{
    public static final String DEFAULT_PYTHON_VERSION = "3.9.4";

    private static final File PYENV_DIR = new File(System.getProperty("user.home") + "/.pyenv");
    private static final File PYENV_PATH = new File(PYENV_DIR + "/libexec/pyenv");

    public final String mPythonVersion;
    public final File mVirtualEnvName;

    public PythonEnv(String pythonVersion, String virtualEnvName)
    {
        mVirtualEnvName = new File(virtualEnvName);
        mPythonVersion = pythonVersion;
    }

    public File virtualEnvPath() { return new File(String.format("%s/versions/%s/envs/%s", PYENV_DIR, mPythonVersion, mVirtualEnvName)); }

    public File pythonPath() { return new File(virtualEnvPath() + "/bin/python"); }

    @VisibleForTesting
    void installPyenv()
    {
        if(!PYENV_DIR.exists())
        {
            String command;

            command = "curl https://pyenv.run | bash"; // This also installs pyenv-virtualenv to $HOME/.pyenv/plugins/pyenv-virtualenv
            CUP_LOGGER.info("Installing pyenv with command: " + command);
            new BashCommand(command).logLevel(Level.DEBUG).run();
        }
        else
        {
            CUP_LOGGER.warn("Skipping installing pyenv as it exists at: " + PYENV_DIR);
        }

        if(!PYENV_DIR.exists())
        {
            CUP_LOGGER.error("Failed to install pyenv to " + PYENV_DIR);
            System.exit(1);
        }
    }

    @VisibleForTesting
    File getRcFilePath() throws FileNotFoundException
    {
        String[] basenamesToTry = {".zshrc", ".bashrc"};

        String pathToRcFile = null;
        for(String basename : basenamesToTry)
        {
            String path = System.getProperty("user.home") + "/" + basename;
            if(new File(path).isFile())
                pathToRcFile = path;
        }

        if(pathToRcFile == null)
            throw new FileNotFoundException("Failed to locate one of the following rc files: " + String.join(", ", basenamesToTry));

        return new File(pathToRcFile);
    }

    @VisibleForTesting
    void addPyenvPathsToRcFile()
    {
        try
        {
            File RcFilePath = getRcFilePath();

            String lines = String.join("\n",
                "# Pyenv paths",
                "export PYENV_ROOT="+ PYENV_DIR,
                "[[ -d $PYENV_ROOT/bin ]] && export PATH=\"$PYENV_ROOT/bin:$PATH\"",
                "eval \"$(pyenv init -)\"",
                "eval \"$(pyenv virtualenv-init -)\""
            );

            boolean linesExist = FileUtils.readFileToString(RcFilePath, "UTF-8").contains(lines);
            if(!linesExist)
            {
                CUP_LOGGER.info("Adding pyenv paths to rc file: " + RcFilePath);
                String command = String.format("echo -e '\n%s' >> %s", lines, RcFilePath);
                new BashCommand(command).logLevel(null).run();
            }
            else
            {
                CUP_LOGGER.warn("Skipping adding pyenv paths to rc file: " + RcFilePath);
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.warn("Failed to add pyenv paths to rc file (using pyenv in interactive shell will not work): " + e);
        }
    }

    @VisibleForTesting
    void installPython()
    {
        File pythonPathInPyenv = new File(String.format("%s/versions/%s/bin/python", PYENV_DIR, mPythonVersion));

        if(!pythonPathInPyenv.exists())
        {
            String command = String.format("%s install %s", PYENV_PATH, mPythonVersion);
            CUP_LOGGER.info("Installing python version {} with command: {}", mPythonVersion, command);
            new BashCommand(command).logLevel(Level.DEBUG).run();
        }
        else
        {
            CUP_LOGGER.warn("Skipping installing python {} as it exists at: {}", mPythonVersion, pythonPathInPyenv);
        }
    }

    @VisibleForTesting
    void createVirtualEnvironment()
    {
        if(!pythonPath().exists())
        {
            String command = String.format("%s virtualenv %s %s", PYENV_PATH, mPythonVersion, mVirtualEnvName);
            CUP_LOGGER.info("Creating python virtual environment with command: " + command);
            new BashCommand(command).logLevel(Level.DEBUG).run();
        }
        else
        {
            CUP_LOGGER.warn("Skipping creating python virtual environment as binary exists at: " + pythonPath());
        }
    }

    @VisibleForTesting
    void removeExistingPyenv()
    {
        CUP_LOGGER.info("Removing existing pyenv at: " + PYENV_DIR);
        try {
            FileUtils.deleteDirectory(PYENV_DIR);
        } catch(IOException e) {
            CUP_LOGGER.warn("Failed to remove pyenv: " + e);
        }
    }

    @VisibleForTesting
    PythonEnv initialize(boolean removePyenv)
    {
        if(removePyenv)
            removeExistingPyenv();

        if(pythonPath().exists())
            return this;

        installPyenv();
        addPyenvPathsToRcFile();
        installPython();
        createVirtualEnvironment();
        return this;
    }

    public PythonEnv initialize()
    {
        return initialize(false);
    }

    private String pipList()
    {
        return new PythonEnvCommand(this, "pip --disable-pip-version-check list")
                .logLevel(null).run().getStdoutAsString();
    }

    public boolean packageInstalled(String packageName)
    {
        String stdout = pipList();
        return stdout.contains(packageName);
    }

    public boolean packagesInstalled(String[] packageNames)
    {
        String stdout = pipList();

        List<String> missingPackages = new ArrayList<>();
        for(String packageName : packageNames)
        {
            if(!stdout.contains(packageName))
                missingPackages.add(packageName);
        }

        if(missingPackages.size() > 0)
        {
            CUP_LOGGER.warn("Python environment is missing the following packages: " + String.join(", ", missingPackages));
            return false;
        }

        return true;
    }

    public void pipUpgrade()
    {
        String command = "pip install --upgrade pip";
        CUP_LOGGER.info("Upgrading pip with command: " + command);
        new PythonEnvCommand(this, command).showCommand().logLevel(Level.DEBUG).run();
    }

    public void pipInstall(String args)
    {
        String command = "pip install " + args;
        CUP_LOGGER.info("Installing packages with command: " + command);
        new PythonEnvCommand(this, command).showCommand().logLevel(Level.DEBUG).run();
    }
}