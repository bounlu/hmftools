package com.hartwig.hmftools.orange.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeConfig;
import com.hartwig.hmftools.datamodel.orange.OrangeConfig;
import com.hartwig.hmftools.orange.TestOrangeConfigFactory;

import org.junit.Test;

public class OrangeAlgoTest {

    @Test
    public void canRunReportFromTestDirPanel() throws IOException {
        OrangeConfig config = TestOrangeConfigFactory.createPanelConfig();
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canRunReportFromTestDirWGSTumorOnly() throws IOException {
        OrangeConfig config = TestOrangeConfigFactory.createWGSConfigTumorOnly();
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canRunReportFromTestDirWGSTumorNormal() throws IOException {
        OrangeConfig config = TestOrangeConfigFactory.createWGSConfigTumorNormal();
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canRunReportFromTestDirWGTSTumorNormal() throws IOException {
        OrangeConfig config = TestOrangeConfigFactory.createWGTSConfigTumorNormal();
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canCreateReportWithoutTumorDoids() throws IOException {
        OrangeConfig config = ImmutableOrangeConfig.builder()
                .from(TestOrangeConfigFactory.createWGSConfigTumorNormal())
                .primaryTumorDoids(Sets.newHashSet())
                .build();

        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }
}