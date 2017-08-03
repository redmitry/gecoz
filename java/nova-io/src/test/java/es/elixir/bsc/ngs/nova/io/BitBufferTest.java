package es.elixir.bsc.ngs.nova.io;

import java.io.IOException;
import java.nio.ByteBuffer;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Dmitry Repchevsky
 */
public class BitBufferTest {
    
    @Test
    public void writeTest() throws IOException {
        ByteBuffer bb = ByteBuffer.allocate(128);
        BitBuffer buf = new BitBuffer(bb);
        
        buf.writeBits(0b00100011110010101010, 18);
        buf.flush();
        
        Assert.assertEquals(0b10101010, bb.get(0) & 0xFF);
        Assert.assertEquals(0b00111100, bb.get(1) & 0xFF);
        Assert.assertEquals(0b00000010, bb.get(2) & 0xFF);
    }
    
    @Test
    public void readTest() throws IOException {
        ByteBuffer bb = ByteBuffer.allocate(128);
        BitBuffer buf = new BitBuffer(bb);
        
        buf.writeBits(0b00100011110010101010, 18);
        
        buf.rewind();
        
        Assert.assertEquals(0b1010, buf.readBits(4) & 0b1111);
    }

    @Test
    public void flushTest() throws IOException {
        ByteBuffer bb = ByteBuffer.allocate(3);
        BitBuffer buf = new BitBuffer(bb);
        
        buf.writeBits(0b100011110010101010, 18);
        
        buf.flush();
        
        Assert.assertEquals(0b10101010, bb.get(0) & 0xFF);
        Assert.assertEquals(0b00111100, bb.get(1) & 0xFF);
        Assert.assertEquals(0b10, bb.get(2) & 0xFF);
    }
    
    @Test
    public void alignedFlushTest() throws IOException {
        ByteBuffer bb = ByteBuffer.allocate(3);
        BitBuffer buf = new BitBuffer(bb);
        
        buf.writeBits(0b100011110010101010, 18);
        
        buf.flush();
        
        Assert.assertEquals(0b10101010, bb.get(0) & 0xFF);
        Assert.assertEquals(0b00111100, bb.get(1) & 0xFF);
        Assert.assertEquals(0b10, bb.get(2) & 0xFF);
    }
}
