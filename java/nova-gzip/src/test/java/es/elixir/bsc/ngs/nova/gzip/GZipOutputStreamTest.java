package es.elixir.bsc.ngs.nova.gzip;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.Test;

/**
 * @author Dmitry Repchevsky
 */

public class GZipOutputStreamTest {
    
    @Test
    public void test_01() throws IOException {

        byte[] file = Files.readAllBytes(Paths.get("D:/WORK/GENCOM/TEESSTTTT/english50/english.50MB"));
        
        final long start = System.nanoTime();
        try (GZipOutputStream gzip = new GZipOutputStream(new BufferedOutputStream(
                Files.newOutputStream(Paths.get("D:/WORK/GENCOM/TEESSTTTT/english50/english50.test.gzip"),
                StandardOpenOption.WRITE, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)))) {

//                for (int i = 0; i < file.length; i++) {
//                    gzip.write(file[i] & 0xFF);
//                }

                  gzip.write(file);
                
//                int off = 0;
//                int block = 1024 * 32;
//                for (; off < file.length; off += block) {
//                    gzip.write(file, off, Math.min(file.length - off, block));
//                    gzip.flush();
//                }

        } catch (FileNotFoundException ex) {
            Logger.getLogger(GZipOutputStreamTest.class.getName()).log(Level.SEVERE, null, ex);
        }
        final long end = System.nanoTime();
        
        System.out.println("compresssed in " + (end - start)/1000000);
        
        GZipFileInputStream in = new GZipFileInputStream(Paths.get("D:/WORK/GENCOM/TEESSTTTT/english50/english50.test.gzip"));
        
        for (int i = 0; i < file.length; i++) {
            int ch = in.read();
            if ((file[i] & 0xFF) != ch) {
                System.out.println(i + " " + (char)file[i] + "  " + (char)ch);
            }
        }
    }
}
