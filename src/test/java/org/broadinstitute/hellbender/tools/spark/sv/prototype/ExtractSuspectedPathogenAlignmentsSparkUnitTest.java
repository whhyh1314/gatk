package org.broadinstitute.hellbender.tools.spark.sv.prototype;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ExtractSuspectedPathogenAlignmentsSparkUnitTest extends BaseTest {

    @DataProvider(name = "testAlignmentsCoverage")
    private Object[][] CreateTestDataForestAlignmentsCoverage() {
        final Object[][] data = new Object[2][];

        AlignmentInterval one = new AlignmentInterval(new SimpleInterval("chr19", 108240, 108401), 1, 160, TextCigarCodec.decode("160M"), false, 0, 0, 160, false);
        AlignedContig tig = new AlignedContig("asm000000:tig00000", "GAATAATCTATGGCATGAAAGATTTTACCTGTCAACAGTGGCTGGCTCTTCATGGTTGCTACAATGAGTGTGTAAGATTCTGAAGAACTCCTTTAATAAGCCTAAACTTAATGTTCAACTTAGAATAAATACAATTCTTCTAAATTTTTTTGAATAATTT".getBytes(),
                Arrays.asList(one));
        data[0] = new Object[]{tig, 160};

        one = new AlignmentInterval(new SimpleInterval("chr5", 66018509, 66018581), 180, 248, TextCigarCodec.decode("179H35M4D34M7H"), true, 0, 5, 44, false);
        AlignmentInterval two = new AlignmentInterval(new SimpleInterval("chr1", 66379, 66452), 138, 211, TextCigarCodec.decode("137S74M44S"), false, 60, 0, 74, false);
        AlignmentInterval three = new AlignmentInterval(new SimpleInterval("chr1", 66373, 66483), 1, 116, TextCigarCodec.decode("41M5I70M139H"), false, 9, 11, 60, false);
        tig = new AlignedContig("asm000000:tig00012", "TATATATTATTATATAATATAATATATATTATATAATATATTTTATTATATAATATAATATATATTATATAATATAATATATTTTATTATATAAATATATATTATATTATATAATATAATATATATTTATATAATATAAATATATATTATATTATATAATATAATATATATTATATAATATATTTTATTATATAAATATATATTATATTATATAATATAATATATATTTTATTATATAATATATATTATATATTT".getBytes(),
                Arrays.asList(three, two, one));
        data[1] = new Object[]{tig,227};

        return data;
    }

    @Test(dataProvider = "testAlignmentsCoverage", groups = "sv")
    public void testAlignmentsCoverage(final AlignedContig input, final int expectedCoverage) {
        Assert.assertEquals(ExtractSuspectedPathogenAlignmentsSpark.alignmentsCoverage(input), expectedCoverage);
    }

}
