package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public final class ExampleReadWalkerWithVariantsIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/examples/";

    @Test
    public void testExampleReadWalkerWithVariants() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -L 1:100-200" +
                " -R " + hg19MiniReference +
                " -I tumor:" + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                " -V " + TEST_DATA_DIRECTORY + "feature_data_source_test.vcf" +
                " -O %s",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleReadWalkerWithVariantsIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleReadWalkerWithVariants", this);
    }
}