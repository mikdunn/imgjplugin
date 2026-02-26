package fftanalysis.imagej;

import ij.IJ;
import ij.ImagePlus;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.nio.file.*;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.*;

public class GenerateTileMontageFromSourceTest {

    @Test
    void runTileMontageFromConfiguredSourceImage() throws Exception {
        String input = System.getProperty("fiba.input", "C:/Users/dunnmk/Downloads/C15D5P001 (1).jpg");
        String output = System.getProperty("fiba.output");
        String base = System.getProperty("fiba.base", "C15D5P001_1");
        int tilesY = Integer.parseInt(System.getProperty("fiba.tilesY", "10"));

        if (output == null || output.trim().isEmpty()) {
            throw new IllegalArgumentException("Missing -Dfiba.output=<absolute output directory>");
        }

        Path inputPath = Paths.get(input).toAbsolutePath().normalize();
        Path outputDir = Paths.get(output).toAbsolutePath().normalize();

        assertTrue(Files.exists(inputPath), "Source image does not exist: " + inputPath);
        Files.createDirectories(outputDir);

        ImagePlus imp = IJ.openImage(inputPath.toString());
        assertNotNull(imp, "ImageJ could not open source image: " + inputPath);
        assertNotNull(imp.getProcessor(), "Opened image has null processor: " + inputPath);
        imp.setTitle(base + ".jpg");

        FIBA_Tile_Montage plugin = new FIBA_Tile_Montage();
        String opts = "outputDir=" + outputDir
                + " overwrite=true"
                + " tilesY=" + tilesY
                + " distributeYToEdges=true"
                + " trackX=true"
                + " saveOverlay=true"
                + " saveMontage=true"
                + " savePlot=true"
                + " savePerTile=true"
                + " saveSol=true"
                + " saveCsv=true"
                + " alpha=0.3 beta=0.3 gamma=0.3 rmin=4";

        plugin.setup(opts, imp);
        plugin.run(imp.getProcessor());

        Path overlay = outputDir.resolve(base + "_tile_boxes.jpg");
        Path stack = outputDir.resolve(base + "_tile_montage.jpg");
        Path profile = outputDir.resolve(base + "_tile_profile.jpg");
        Path csv = outputDir.resolve(base + "_tile_results.csv");

        assertTrue(Files.isRegularFile(overlay), "Missing generated overlay: " + overlay);
        assertTrue(Files.isRegularFile(stack), "Missing generated tile stack: " + stack);
        assertTrue(Files.isRegularFile(profile), "Missing generated tile profile: " + profile);
        assertTrue(Files.isRegularFile(csv), "Missing generated CSV: " + csv);

        for (int i = 1; i <= tilesY; i++) {
            assertTrue(Files.isRegularFile(outputDir.resolve(base + "_tile" + i + "_crop.jpg")), "Missing crop tile " + i);
            assertTrue(Files.isRegularFile(outputDir.resolve(base + "_tile" + i + "_fft.jpg")), "Missing fft tile " + i);
            assertTrue(Files.isRegularFile(outputDir.resolve(base + "_tile" + i + "_polar.jpg")), "Missing polar tile " + i);
            assertTrue(Files.isRegularFile(outputDir.resolve(base + "_tile" + i + "_sol.jpg")), "Missing SOL tile " + i);
            assertTrue(Files.isRegularFile(outputDir.resolve(base + "_tile" + i + "_mask.jpg")), "Missing mask tile " + i);
            assertTrue(Files.isRegularFile(outputDir.resolve(base + "_tile" + i + "_rec.jpg")), "Missing reconstruction tile " + i);
        }

        try (Stream<Path> s = Files.list(outputDir)) {
            List<String> generated = s
                    .map(p -> p.getFileName().toString())
                    .filter(n -> n.startsWith(base + "_"))
                    .sorted()
                    .collect(Collectors.toList());
            assertFalse(generated.isEmpty(), "No generated files found for base " + base);
            System.out.println("Generated files for " + base + ":");
            for (String n : generated) {
                System.out.println(" - " + n);
            }
        }
    }
}
