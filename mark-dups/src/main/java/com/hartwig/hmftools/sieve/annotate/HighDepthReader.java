package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

// TODO(m_cooper): Code duplication?
public class HighDepthReader
{
    private static final String DELIMITER = ",";
    private static final int NUM_FIELDS = 4;

    static public List<HighDepthRegion> readFromFile(final String filepath)
    {
        final List<HighDepthRegion> records = new ArrayList<>();

        try
        {
            final BufferedReader reader = new BufferedReader(new FileReader(filepath));

            // Drop header.
            reader.readLine();

            String line;
            while((line = reader.readLine()) != null)
            {
                final String[] fields = line.split(DELIMITER);
                if(fields.length != NUM_FIELDS)
                {
                    MD_LOGGER.error("High depth file {} contains a record with {} fields. There should only be {} fields. ", filepath, fields.length, NUM_FIELDS);
                    System.exit(1);
                }

                final HighDepthRegion record = parseFields(fields, filepath);
                records.add(record);
            }

            reader.close();

        }
        catch(Exception e)
        {
            MD_LOGGER.error("An exception was raised while reading the high depth file {}: {}", filepath, e.toString());
            System.exit(1);
        }

        return records;
    }

    static private HighDepthRegion parseFields(final String[] fields, final String filepath)
    {
        String chromosome = fields[0];
        int posStart = parseInt(fields[1], "PosStart", filepath);
        int posEnd = parseInt(fields[2], "PosEnd", filepath);
        int baseDepth = parseInt(fields[3], "BaseDepth", filepath);

        return new HighDepthRegion(chromosome, posStart, posEnd, baseDepth);
    }

    static private int parseInt(final String str, final String fieldName, final String filepath)
    {
        try
        {
            return Integer.parseInt(str);
        }
        catch(NumberFormatException e)
        {
            MD_LOGGER.error("While reading a record in the high depth file {}, failed to parse the '{}' field to an int.", filepath, fieldName);
            System.exit(1);
        }

        return 0;
    }
}