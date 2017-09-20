/**
 * *****************************************************************************
 * Copyright (C) 2015 Spanish National Bioinformatics Institute (INB) and
 * Barcelona Supercomputing Center
 *
 * Modifications to the initial code base are copyright of their respective
 * authors, or their employers as appropriate.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 *****************************************************************************
 */

package es.elixir.bsc.ngs.nova.algo.string;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.InvalidPathException;
import java.nio.file.Path;
import java.util.Arrays;

/**
 * This class implements SACA-K suffix array construction method (doi:10.1145/2493175.2493180).
 * When there is enough space for buckets may be allocated, the iteration is done via SA-IS algorithm (doi:10.1109/DCC.2009.42).
 * 
 * @author Dmitry Repchevsky
 */

public class SuffixArray {
    
    /* the cue-card
     * 
     *  s   - the original array
     *  sa  - the suffix array (where the all sorting is performed)
     *  k   - the starting pointer to the reduced string
     *  m   - the end pointer of the reduced string + 1 ( sa[k .. m - 1] )
     *  n   - the number of LMS suffixes ( lms[0 .. n - 1 )
     *  l   - the alphabet size (or '0' if all substrings are equal)
     *  c   - the counter of symbols in the string (either c[] in level 0, or the pointer to SA if level > 0)
     *  bk1 - the 'L' buckets
     *  bk2 - the 'S' buckets
     */
    public static void main(String[] args) {
        if (args.length == 0) {
            System.out.println("input file missed.\nusage: SuffixArray input [output]");
        }
        
        try {
            Path input = FileSystems.getDefault().getPath(args[0]);
            byte[] s = Files.readAllBytes(input);
            int[] sa = new int[s.length];
            final long start = System.currentTimeMillis();
            SuffixArray.suffix(ByteBuffer.wrap(s), sa);
            final long end = System.currentTimeMillis();
            System.out.println("finished " + s.length + " in " + (end - start) + " ms.");
            
            if (args.length == 2) {
                PrintWriter writer = new PrintWriter(args[1]);
                for (int i = 0; i < sa.length; i++) {
                    writer.write(Integer.toString(sa[i]));
                    writer.write((i + 1) % 10 == 0 ? '\n' : '\t');
                }
                writer.close();
            }
        } catch (InvalidPathException ex) {
            System.out.println("wrong input file path:\n" + args[0]);
        } catch (FileNotFoundException ex) {
            System.out.println("cant create file:\n" + args[1]);
        } catch(IOException ex) {
            System.out.println(ex.getMessage());
        }
    }
    
    public static int[] suffix(ByteBuffer s) {
        int[] sa = new int[s.limit()];
        suffix(s, sa);
        return sa;
    }

    /**
     * 
     * @param s the original array
     * @param sa the suffix array
     * 
     * @return alphabet buckets' indexes ( 0 .. 255 ) or null if s.length != sa.length
     */
    public static int[] suffix(ByteBuffer s, int[] sa) {
        if (s.limit() != sa.length) {
            return null;
        }

        int[] c = new int[256];
        count(s, c);
        
        int[] bk1 = new int[256];
        int[] bk2 = new int[256];
        
        buckets(c, bk1, bk2);

        final int n = fillSA(s, sa, Arrays.copyOf(bk2, bk2.length));
        if (n > 1) {
            sortLMS(s, sa, Arrays.copyOf(bk1, bk1.length), Arrays.copyOf(bk2, bk2.length));
            compactLMS(sa, sa.length, n);

            int l = nameSubstrS(s, sa, n);
            if (l > 1) {
                Arrays.fill(sa, 0, n, 0); // clear buckets
                if (l == sa.length - n) {
                    sortSA0_saka_k(sa, sa.length - n, sa.length);
                } else {
                    sortSA0_sais(sa, sa.length - n, sa.length, l);
                }
                getBackLMS(s, sa, n);
            }

            putSuffix(s, sa, n, Arrays.copyOf(bk2, bk2.length));
        }

        induceLMS(s, sa, bk1, Arrays.copyOf(bk2, bk2.length));
        return bk2;
    }

    /**
     * Counts all alphabet symbols in original string 's' to the c[] array.
     * 
     * @param s the original array
     * @param c the counter array to be calculated
     */
    private static void count(final ByteBuffer s, final int[] c) {
        for (int i = 0, n = s.limit(); i < n; i++) {
            c[s.get(i) & 0xFF]++;
        }
    }
    
    /**
     * Calculates starts and ends of the buckets
     * 
     * @param c - number elements in the buckets
     * @param bk1 start pointers of buckets (calculated by the method)
     * @param bk2 end pointers of buckets (calculated by the method)
     */
    private static void buckets(final int[] c, final int[] bk1, final int[] bk2) {
        for (int i = 0, idx = 0; i < c.length; i++) {
            final int p = c[i];
            if (p > 0) {
                bk1[i] = idx;
                idx += p;
                bk2[i] = idx - 1;
            }
        }
    }

    /**
     * Recursive sort of SA0 using inplace bucket counters (saka-k)
     * 
     * @param sa the suffix array
     * @param k the pointer to the first name in the SA0 array
     * @param m the pointer to the last name in the SA0 array (SA0.length)
     */
    private static void sortSA0_saka_k(final int[] sa, final int k, final int m) {
        final int n = fillSA0_saka_k(sa, k, m);
        if (n > 1) {
            sortLMS_saka_k(sa, k, m - k, false);
            compactLMS_saka_k(sa, m - k, n);

            int l = nameSubstrSA(sa, k, m, n);
            if (l == 0) {
                l = k - n;
                for(int i = 0; i < n; i++) {
                    sa[sa[i + l]] = i;
                }            
            } else if (l == 1) {
                induceLMS(sa, k - n, k, n, n, n, l);
            } 
            else if (l == k - n) {
                sortSA0_saka_k(sa, k - n, k);
            } else {
                sortSA0_sais(sa, k - n, k, l);
            }

            getBackLMS(sa, k, m, n);
            putSuffix_saka_k(sa, k, m, n);
        }
        sortLMS_saka_k(sa, k, m - k, true);
    }

    /**
     * Recursive sort of SA0 using a bucket counter array.
     * 
     * @param k the pointer to the first name in the SA0 array
     * @param m the pointer to the last name in the SA0 array (SA0.length)
     * @param c the pointer to the first index of bucket counter array
     */
    private static void sortSA0_sais(final int[] sa, final int k, final int m, final int alph) {
        
        final int c = m - k; // counters array is placed right after the LMS
        final int bk1 = c + alph * 2 < k ? c + alph : c;
        final int bk2 = bk1 + alph * 2 < k ? bk1 + alph : bk1;
        
        final int n = fillSA0_sais(sa, k, m, c, bk1, bk2, alph);

        if (n > 1) {
            sortLMS(sa, k, m, c, bk1, bk2, alph);
            compactLMS(sa, m - k, n);

            Arrays.fill(sa, c, bk2 + alph, 0); // clear buckets

            int l = nameSubstrSA(sa, k, m, n);
            if (l == 0) {
                l = k - n;
                for(int i = 0; i < n; i++) {
                    sa[sa[i + l]] = i;
                }
            } else if (l == 1) {
                induceLMS(sa, k - n, k, n, n, n, l);
            } else if (l == k - n) {
                sortSA0_saka_k(sa, k - n, k);            
            } else {
                sortSA0_sais(sa, k - n, k, l);            
            }
            
            getBackLMS(sa, k, m, n);
            putSuffix(sa, k, m, n, bk1, bk2, alph);
        }
        induceLMS(sa, k, m, c, bk1, bk2, alph);
    }

    private static void getBackLMS(final ByteBuffer s, final int[] sa, final int n) {
        
        int ptr = sa.length;
        boolean t = false;
        for (int i = s.limit() - 2, b0 = s.get(s.limit() - 1) & 0xFF; i >= 0; i--) {
            final int b1 = s.get(i) & 0xFF;
            if (b1 < b0) {
                t = true;
            } else if (b1 > b0 && t) {
                sa[--ptr] = i + 1;
                t = false;
            }
            b0 = b1;
        }

        for (int i = 0; i < n; i++) {
            sa[i] = sa[sa[i] + ptr];
        }
    }

    /**
     * @param sa the suffix array
     * @param k the pointer to the first name in the SA0 array
     * @param m the end of the SA0 array
     * @param n the number of sorted LMS suffixes
     */
    private static void getBackLMS(final int[] sa, final int k, int m, final int n) {
        
        int ptr = k;
        boolean t = false;
        for (int i = m - 2, b0 = sa[m - 1]; i >= k; i--) {
            final int b1 = sa[i];
            if (b1 < b0) {
                t = true;
            } else if (b1 > b0 && t) {
                sa[--ptr] = i - k + 1;
                t = false;
            }
            b0 = b1;
        }
        
        for (int i = 0; i < n; i++) {
            sa[i] = sa[sa[i] + ptr];
        }
    }
    
    /**
     * 
     * @param s
     * @param sa
     * @param n
     * @param bk2 
     */
    private static void putSuffix(final ByteBuffer s, final int[] sa, final int n, int[] bk2) {

        Arrays.fill(sa, n, sa.length, 0);
        int ch0 = s.get(sa[n - 1]) & 0xFF;
        int b = bk2[ch0];
        
        for(int i = n - 1; i >= 0; i--) {   
            final int l = sa[i];
            sa[i] = 0;
            final int ch = s.get(l) & 0xFF;
            if (ch != ch0) {
                bk2[ch0] = b;
                b = bk2[ch0 = ch];
            }
            sa[b--] = l;
        }
    }
    
    /**
     * 
     * @param sa
     * @param k
     * @param m
     * @param n
     * @param bk1
     * @param bk2
     * @param l 
     */
    private static void putSuffix(final int[] sa, final int k, final int m, final int n, final int bk1, final int bk2, int l) { 
        
        buckets(sa, k, m, 0, bk1, bk2, l, true);
        
        Arrays.fill(sa, n, m - k, 0);
        
        int ch0 = sa[k + sa[n - 1]];
        int b = sa[bk2 + ch0];
        for(int i = n - 1; i >= 0; i--) {
            final int p = sa[i]; 
            sa[i] = 0;
            final int ch = sa[k + p];

            if (ch != ch0) {
                sa[bk2 + ch0] = b;
                b = sa[bk2 + ch];
                ch0 = ch;
            }
            sa[b--] = p;
        }
    }

    /**
     * 
     * @param sa
     * @param k the pointer to the first name in the SA0 array
     * @param m the end of the SA0 array
     * @param n the number of sorted LMS suffixes
     * 
     */
    private static void putSuffix_saka_k(final int[] sa, final int k, final int m, final int n) { 
        Arrays.fill(sa, n, m - k, Integer.MIN_VALUE);
        
        for(int i = n - 1, ptr = Integer.MIN_VALUE, last = Integer.MIN_VALUE; i >= 0; i--) {
          int p = sa[i]; 
          sa[i] = Integer.MIN_VALUE;
          int pos = sa[p + k];
          if (pos != last) {
              last = pos;
              ptr = pos;
          }
          sa[ptr--] = p;
        }
    }

    /**
     * 
     * @param sa the suffix array
     * @param k the pointer to the first name in the SA0 array
     * @param m the pointer to the end of the SA0 array (last element + 1)
     * @param c the pointer to alphabet symbols counters
     * @param bk1 the pointer to the 'L'-type buckets' counter array
     * @param bk2 the pointer to the 'S'-type buckets' counter array
     * @param l the size of alphabet
     * @param end calculate buckets' ends if 'true' ('S' buckets), starts otherwise ('L' buckets).
     */
    private static void buckets(final int[] sa, final int k, final int m, int c, int bk1, int bk2, final int l, final boolean end) {
        
        if (c == 0) {
            // counters need to be recalculated
            c = m - k;
            Arrays.fill(sa, c, c + l, 0);
            for (int i = m - 1; i >= k; i--) {
                sa[c + sa[i]]++;
            }
        }

        if (bk1 != bk2) {
            // calculate both 'L' and 'S' buckets' pointers
            for (int idx = 0, n = bk2 + l; bk2 < n; c++, bk1++, bk2++) {
                final int p = sa[c];
                sa[bk1] = idx;
                idx += p;
                sa[bk2] = idx - 1;
            }
        } else if (end) {
            for (int idx = 0, n = bk2 + l; bk2 < n; c++, bk2++) {
                idx += sa[c];
                sa[bk2] = idx - 1;
            }            
        } else {
            for (int idx = 0, n = bk1 + l; bk1 < n; c++, bk1++) {
                final int p = sa[c];
                sa[bk1] = idx;
                idx += p;
            }            
        }
    }
    
    /**
     * Puts all LMS-suffixes into the corresponding buckets
     * 
     * @param s
     * @param sa the suffix array
     * @param bk2 'S' type buckets (right-most index)
     * 
     * @return the number of LMS suffixes
     */
    private static int fillSA(final ByteBuffer s, final int sa[], final int[] bk2) {
        int m = 0;
        
        boolean t = false;
        for (int i = s.limit() - 2, b0 = s.get(s.limit() - 1) & 0xFF; i >= 0; i--) {
            final int b1 = s.get(i) & 0xFF;
            if (b1 < b0) {
                t = true;
            } else if (t && b1 > b0) {
                m++;
                final int ptr = bk2[b0]--;
                sa[ptr] = i + 1;
                t = false;                    
            }
            b0 = b1;
        }
        
        return m;
    }
    
    /**
     * Fills LMS suffixes
     * 
     * @param sa the suffix array
     * @param k
     * @param m
     * @param c
     * @param bk1 the pointer to buckets' start
     * @param bk2 the pointer to buckets' end
     * @param l the alphabet size
     */
    private static int fillSA0_sais(final int sa[], final int k, final int m, final int c, final int bk1, final int bk2, final int l) {
        buckets(sa, k, m, c, bk1, bk2, l, true); // get buckets' ends

        int s = 0;
        boolean t = false;
        for (int i = m - 2, b0 = sa[m - 1]; i >= k; i--) {
            final int b1 = sa[i];
            if (b1 < b0) {
                t = true;
            } else if (b1 > b0 && t) {
                s++;
                sa[sa[sa[i + 1] + bk2]--] = i - k + 1;
                t = false;
            }
            b0 = b1;
        }
        return s;
    }
    
    /**
     * @param sa the suffix array
     * @param k the pointer to the first name in the SA0 array
     * @param m the end of the SA0 array ( SA0.length &lt; SA.length)
     */
    private static int fillSA0_saka_k(final int sa[], int k, int m) {
        Arrays.fill(sa, 0, m - k, Integer.MIN_VALUE); // clean
        
        int s = 0;
        boolean t = false;
        for (int i = m - 2, b0 = sa[m - 1]; i >= k; i--) {
            final int b1 = sa[i];
            
            if (b1 < b0) {
                t = true;
            } else if (b1 > b0 && t) {
                t = false;
                s++;
                int c = sa[b0];
                if (c >= 0) {
                    int ptr = b0;
                    do {
                        int tmp = sa[++ptr];
                        sa[ptr] = c;
                        c = tmp;                            
                    } while (c >= 0 || c == Integer.MIN_VALUE);
                    c = Integer.MIN_VALUE;
                }
                if(c == Integer.MIN_VALUE) {
                    int ptr = b0;
                    if (ptr > 0 && sa[ptr - 1] == Integer.MIN_VALUE) {
                        sa[ptr] = ~--ptr;
                    }
                    sa[ptr] = i - k + 1;
                } else {
                    c = ~c;
                    if (c > 0 && sa[c - 1] == Integer.MIN_VALUE) {
                        sa[b0] = ~--c;
                    } else {
                        System.arraycopy(sa, c, sa, c + 1, b0 - c);
                    }
                    sa[c] = i - k + 1;
                }
            }
            b0 = b1;
        }

        for (int i = m - k; i >= 0; i--) {
            final int l = sa[i];
            if (l < 0 && l != Integer.MIN_VALUE) {
                while (i > ~l) {
                    sa[i] = sa[--i];
                }
                sa[i] = Integer.MIN_VALUE;
            }
        }
        
        return s;
    }
    
    /**
     * Sorts LMS suffixes of the original 's' array
     * 
     * @param s
     * @param sa the suffix array
     * @param bk1 start pointers of the buckets
     * @param bk2 end pointers of the buckets
     * 
     * @return the last index of sorted LMS suffixes
     */
    private static void sortLMS(final ByteBuffer s, final int[] sa, final int[] bk1, final int[] bk2) {
          
        // sort LMS from left to right

        final int l0 = sa.length - 1;
        int ch0 = s.get(l0) & 0xFF;
        int b = bk1[ch0];
        sa[b++] = l0 == 0 || ch0 > (s.get(l0 - 1) & 0xFF) ? ~l0 : l0;
        
        for (int i = 0; i <= l0; i++) {
            int l = sa[i];
            if (l > 0) {
                final int ch = s.get(--l) & 0xFF;
                if (ch != ch0) {
                    bk1[ch0] = b;
                    b = bk1[ch0 = ch];
                }
                sa[b++] = l == 0 || ch > (s.get(l - 1) & 0xFF) ? ~l : l;
                sa[i] = 0;
            } else if (l != 0) {
                sa[i] = ~l;
            }           
        }
        
        // sort LMS from rigth to left

        ch0 = 0;
        b = bk2[0];
        for (int i = l0; i >= 0; i--) {
            final int l = sa[i] - 1;
            if(l > 0) {
                final int ch = s.get(l) & 0xFF;
                if (ch != ch0) {
                    bk2[ch0] = b;
                    b = bk2[ch0 = ch];
                }
                sa[b--] = ch < (s.get(l - 1) & 0xFF) ? ~l : l;
            }
        }
    }

    /**
     * 
     * @param s
     * @param sa
     * @param bk1
     * @param bk2 
     */
    private static void induceLMS(final ByteBuffer s, final int[] sa, final int[] bk1, final int[] bk2) {
        
        final int l0 = sa.length - 1;
        int ch0 = s.get(l0) & 0xFF;
        int b = bk1[ch0];
        sa[b++] = l0 == 0 || ch0 > (s.get(l0 - 1) & 0xFF) ? ~l0 : l0;
        
        for (int i = 0; i <= l0; i++) {
            int p = sa[i];
            if (p > 0) {
                final int ch = s.get(--p) & 0xFF;
                if (ch != ch0) {
                    bk1[ch0] = b;
                    b = bk1[ch0 = ch];
                }
                sa[b++] = p == 0 || ch > (s.get(p - 1) & 0xFF) ? ~p : p;
            } else if (p != 0) {
                sa[i] = ~p;
            }           
        }
        
        // sort LMS from rigth to left

        ch0 = 0;
        b = bk2[0];
        for (int i = l0; i >= 0; i--) {
            int p = sa[i];
            if(p > 0) {
                final int ch = s.get(--p) & 0xFF;
                if (ch != ch0) {
                    bk2[ch0] = b;
                    b = bk2[ch0 = ch];
                }
                if (b < i) {
                    sa[b--] = p == 0 || ch < (s.get(p - 1) & 0xFF) ? ~p : p;
                }
            } else if (p != 0) {
                sa[i] = ~p;
            }
        }
    }
    
    /**
     * 
     * @param sa
     * @param k the pointer to the first name in the SA0 array
     * @param m the pointer to the end of the SA0 array
     * @param bk1 the pointer to the bucket array
     * @param suffix 
     */
    private static void sortLMS(final int[] sa, final int k, final int m, final int c, final int bk1, final int bk2, final int alph) {
        
        if (c == bk1) {
            buckets(sa, k, m, 0, bk1, bk1, alph, false);
        } else if (bk1 == bk2) {
            buckets(sa, k, m, c, bk1, bk1, alph, false);
        }
        
        final int n = m - k;
        
        // sort LMS from left to right

        int p = m - 1;
        int ch0 = sa[p];
        int b = sa[bk1 + ch0];
        sa[b++] = p == 0 || ch0 > sa[p - 1] ? ~(p - k) : p - k;

        for (int i = 0; i < n; i++) {
            p = sa[i];
            if (p > 0) {
                final int ch = sa[--p + k];
                if (ch != ch0) {
                    sa[bk1 + ch0] = b;
                    b = sa[bk1 + ch];
                    ch0 = ch;
                }
                sa[b++] = p == 0 || ch > sa[k + p - 1] ? ~p : p;
                sa[i] = 0;
            } else if (p != 0) {
                sa[i] = ~p;
            }
        }
        

        buckets(sa, k, m, c == bk1 ? 0 : c, bk2, bk2, alph, true); // get buckets' ends
        
        // sort LMS from rigth to left

        ch0 = 0;
        b = sa[bk2];

        for (int i = n - 1; i >= 0; i--) {
            p = sa[i] - 1;
            if(p > 0) {
                final int ch = sa[k + p];
                
                if (ch != ch0) {
                    sa[bk2 + ch0] = b;
                    b = sa[bk2 + ch];
                    ch0 = ch;
                }
                
                sa[b--] = ch < sa[k + p - 1] ? ~p : p;
            }
        }
    }

    private static void induceLMS(final int[] sa, final int k, final int m, final int c, final int bk1, final int bk2, final int alph) {
        
        buckets(sa, k, m, c == bk1 ? 0 : c, bk2, bk2, alph, bk1 != bk2);
        
        final int n = m - k;
        
        // sort LMS from left to right

        int p = m - 1;
        int ch0 = sa[p];
        int b = sa[bk1 + ch0];
        sa[b++] = p == 0 || ch0 > sa[p - 1] ? ~(p - k) : p - k;
        
        for (int i = 0; i < n; i++) {
            p = sa[i];
            if (p > 0) {
                final int ch = sa[--p + k];
                if (ch != ch0) {
                    sa[bk1 + ch0] = b;
                    b = sa[bk1 + ch];
                    ch0 = ch;
                }
                sa[b++] = p == 0 || ch > sa[k + p - 1] ? ~p : p;
            } 
            else if (p != 0) {
                sa[i] = ~p;
            }           
        }
        
        if (bk1 == bk2) {
            // there was no space for bk2
            // probably there was no space for bk1 so c must be also recalculated (pass 0)
            buckets(sa, k, m, n == bk1 ? 0 : n, bk1, bk2, alph, true); // get buckets' ends
        }
        
        // sort LMS from rigth to left
        
        ch0 = 0;
        b = sa[bk2];
        
        for (int i = m - k - 1; i >= 0; i--) {
            p = sa[i];
            if(p > 0) {
                final int ch = sa[--p + k];
                
                if (ch != ch0) {
                    sa[bk2 + ch0] = b;
                    b = sa[bk2 + ch];
                    ch0 = ch;
                }

                if (i > b) {
                    sa[b--] = p == 0 || ch < sa[k + p - 1] ? ~p : p;
                }
            } else if (p != 0) {
                sa[i] = ~p;
            }
        }
    }

    /**
     * Compacts LMS suffixes left.
     * 
     * @param sa the suffix array
     * @param m the SA0 length
     * @param n the number of LMS suffixes
     * 
     * @return the number of compacted 'S*' suffixes (the last suffix + 1)
     */
    private static void compactLMS(final int[] sa, final int m, final int n) {
        int idx = 0;
        for (int i = 0; i < n; i++) {
            final int p = sa[i];
            if (p < 0) {
                sa[idx++] = ~p;
            }
        }
        for (int i = n; i < m; i++) {
            final int p = sa[i];
            sa[i] = 0;
            if (p < 0) {
                sa[idx++] = ~p;
            }
        }
    }

    private static void compactLMS_saka_k(final int[] sa, final int m, final int n) {
        int idx = 0;
        for (int i = 0; i < n; i++) {
            final int p = sa[i];
            if (p > 0) {
                sa[idx++] = p;
            }
        }
        for (int i = n; i < m; i++) {
            final int p = sa[i];
            sa[i] = 0;
            if (p > 0) {
                sa[idx++] = p;
            }
        }
    }
    
    /**
     * Sorts LMS suffixes in the SA0 array using SAKA-K algorithm (using inplace buckets' counters).
     * 
     * @param sa the suffix array
     * @param k the pointer to the first element in the SA0 array
     * @param n the number of LMS suffixes
     * @param suffix whether to sort or create suffix array
     */
    private static void sortLMS_saka_k(final int[] sa, final int k, int n, final boolean suffix) {
        
        // sort LMS from left to right

        n--;
        int ptr = sa[n + k];
        if (ptr < n && sa[ptr + 1] == Integer.MIN_VALUE) {
            sa[ptr] = ~++ptr;
        }
        sa[ptr] = n;
                        
        for (int i = 0; i <= n; i++) {
            final int p = sa[i];
            if (p > 0) {
                ptr = sa[p + k - 1];
                final int ptr2 = sa[p + k];
                
                if (ptr >= ptr2) {

                    int c = sa[ptr];                    
                    if (c >= 0) {
                        // occupied by a previous bucket
                        int _ptr = ptr;
                        do {
                            int tmp = sa[--_ptr];
                            sa[_ptr] = c;
                            c = tmp;                            
                        } while (c >= 0 || c == Integer.MIN_VALUE); // untill the previous bucket counter found
                        if (i > _ptr) { // && i <= ptr) {
                            i--;
                        }
                        c = Integer.MIN_VALUE;
                    }
                    if (c == Integer.MIN_VALUE) {
                        if (ptr < n && sa[ptr + 1] == Integer.MIN_VALUE) {
                            sa[ptr] = ~++ptr;
                        }
                        sa[ptr] = p - 1;
                    } else {
                        c = ~c; // is a bucket pointer
                        if (c < n && sa[c + 1] == Integer.MIN_VALUE) {
                            sa[ptr] = ~++c;
                        } else {
                            System.arraycopy(sa, ptr + 1, sa, ptr, c - ptr);

                            if (i > ptr) {// && i <= c) {
                              i--;
                            }
                        }
                        sa[c] = p - 1;
                    }
                    if (suffix) {
                        // clean if 'S*'
                        if (p < n) {
                            final int ptr3 = sa[p + k + 1];
                            if (ptr2 < ptr3 || (ptr2 == ptr3 && ptr2 >= i)) {
                                sa[i] = Integer.MIN_VALUE;
                            }
                        }
                    } else {
                        sa[i] = Integer.MIN_VALUE;
                    }
                }
            }            
        }
        
        for (int i = 0; i <= n; i++) { // i = 0
            final int p = sa[i];
            if (p < 0 && p != Integer.MIN_VALUE) {
                while (i < ~p) {
                    sa[i] = sa[++i];
                }
                sa[i] = Integer.MIN_VALUE;
            }
        }
        
        // sort LMS from rigth to left

        for (int i = n; i >= 0; i--) {
            final int p = sa[i];
            if (p > 0) {
                ptr = sa[p + k - 1];
                final int ptr2 = sa[p + k];
                if (ptr < ptr2 || (ptr == ptr2 && ptr > i)) {
                    int c = sa[ptr];
                    if (c >= 0) {
                        int _ptr = ptr;
                        do {
                            int tmp = sa[++_ptr];
                            sa[_ptr] = c;
                            c = tmp;                            
                        } while (c >= 0 || c == Integer.MIN_VALUE);
                        if (i < _ptr) { // && i >= ptr) {
                            i++;
                        }
                        c = Integer.MIN_VALUE;
                    }
                    if(c == Integer.MIN_VALUE) {
                        if (ptr > 0 && sa[ptr - 1] == Integer.MIN_VALUE) {
                            sa[ptr] = ~--ptr;
                        }
                        sa[ptr] = p - 1;
                    } else {
                        c = ~c;
                        if (c > 0 && sa[c - 1] == Integer.MIN_VALUE) {
                            sa[ptr] = ~--c;
                        } else {
                            System.arraycopy(sa, c, sa, c + 1, ptr - c);
                            if (i < ptr) { // && i >= c) {
                              i++;
                            }
                        }
                        sa[c] = p - 1;
                    }
                    if (!suffix) {
                       sa[i] = Integer.MIN_VALUE; 
                    }
                }
            }
        }

        if (!suffix) {
            for (int i = n; i >= 0; i--) {
                final int p = sa[i];
                if (p < 0 && p != Integer.MIN_VALUE) {
                    while (i > ~p) {
                        sa[i] = sa[--i];
                    }
                    sa[i] = Integer.MIN_VALUE;
                }
            }
        }
    }

    /**
     * 
     * @param n the last index of LMS suffixes in the SA (0 .. n)
     * 
     * @return the first LMS names index in the SA (index .. SA.length)
     */
    private static int nameSubstrS(final ByteBuffer s, final int[] sa, final int n) {
        
        if (n * 3 < sa.length) {
            return nameSubstrS_sais(s, sa, n);
        }
        
        final int l = nameSubstrS_saka_k(s, sa, n);

        return l <= 1 ? l : compactNames_saka_k(sa, sa.length);
    }

    private static int nameSubstrS_sais(final ByteBuffer s, final int[] sa, final int n) {
        
        final int h = sa.length >> 1;
        int alph = -1;

        for (int i = 0, ptr1 = 0, len1 = 0; i < n; i++) {
            final int ptr2 = sa[i];
            final int pos = h + (ptr2 >> 1);
            int len2 = getNextLMS(s, ptr2) - ptr2 + 1;
            if (len1 != len2) {
                len1 = len2;
            } else {
                for (int j = ptr1, k = ptr2; len2 > 0; j++, k++, len2--) {
                    if (s.get(j) != s.get(k)) {
                        break;
                    }
                }
            }

            sa[pos] = len2 == 0 ? ~alph : ~++alph;

            ptr1 = ptr2;
        }
        
        if (alph++ == 0 || alph == n) {
            // either all equals (alphabet = 0) or different
            Arrays.fill(sa, h, sa.length, 0); // clean
            return 0;
        }

        int idx = sa.length;
        if (n + n <= h) {
            // count and compact names at the same time
            for (int i = sa.length - 1; i >= h; i--) {
                int ptr = sa[i];
                if (ptr < 0) {
                    sa[i] = 0;
                    ptr = ~ptr;
                    sa[--idx] = ptr;
                    sa[n + ptr]++; // counter
                }
            }
        } else {
            // compact names
            for (int i = sa.length - 1; i >= h; i--) {
                final int ptr = sa[i];
                if (ptr < 0) {
                    sa[i] = 0;
                    sa[--idx] = ~ptr;
                }
            }
            // count names
            for (int i = idx; i < sa.length; i++) {
                sa[n + sa[i]]++;
            }
        }

        return alph;
    }
    
    /**
     * <p>
     * Naming LMS substring according SACA-K.
     * <br/>
     * If LSM suffixes are equal or all different, no changes performed.
     * </p>
     * 
     * @param s the original string to calculate the suffix array for.
     * @param sa the array that contains sorted LMS suffixes.
     * @param n the number of LMS suffixes in the sa (0 .. n - 1).
     * 
     * @return 
     */
    private static int nameSubstrS_saka_k(final ByteBuffer s, final int[] sa, final int n) {
        
        final int h = sa.length >> 1; // part by half
        int alph = -1;

        for (int i = 0, ptr1 = 0, len1 = 0, bk = 0; i < n; i++) {
            final int ptr2 = sa[i];
            final int pos = h + (ptr2 >> 1);

            int len2 = getNextLMS(s, ptr2) - ptr2 + 1;
            if (len1 != len2) {
                len1 = len2;
            } else {
                for (int j = ptr1, p = ptr2; len2 > 0; j++, p++, len2--) {
                    if (s.get(j) != s.get(p)) {
                        break;
                    }
                }
            }

            if (len2 == 0) { // equal substrings
                sa[pos] = ~bk;
            } else {
                sa[pos] = ~i;
                if (i > bk + 1) {
                    sa[bk] = bk - i;
                    sa[++bk] = ~alph;
                }
                alph++;
                bk = i;
            }

            ptr1 = ptr2;
        }

        return ++alph == n ? 0 : alph;
    }

    /**
     * 
     * @param k the pointer to the first name in the SA0 array
     * @param m the pointer to the end of the SA0 array
     * @param n the number of sorted LMS suffixes
     * 
     * @return - LMS names index in the SA0 (index .. n).
     *           Returns 0 if all names are different.
     */
    private static int nameSubstrSA(final int[] sa, int k, final int m, final int n) {

        if (n * 3 < k) {
            return nameSubstrSA_sais(sa, k, m, n);
        }
        
        final int l = nameSubstrSA_saka_k(sa, k, m, n);
        
        if (n <= k - n - l) {
            // compact LMS names and set them usual way (0 .. alph)
            compactNames_sais(sa, k, n, l);
            Arrays.fill(sa, 0, n, 0); // compactNames_sais() leaves garbage
            return l == n ? 0 : l;
        }
        
        // compact LMS names and set up 'S' buckets
        final int idx = compactNames_saka_k(sa, k);
        return l == n ? 0 : idx;
    }

    private static int nameSubstrSA_sais(final int[] sa, int k, final int m, final int n) {

        final int h = k >> 1; // part by half
        
        int alph = -1;
        
        for (int i = 0, ptr1 = 0, len1 = 0; i < n; i++) {
            final int ptr2 = sa[i];
            sa[i] = 0;
            final int pos = h + (ptr2 >> 1);
            int len2 = getNextLMS(sa, ptr2 + k, m - 1) - ptr2 - k + 1;
            if (len1 != len2) {
                len1 = len2;
            } else {
                for (int j = ptr1 + k, p = ptr2 + k; len2 > 0; j++, p++, len2--) {
                    if (sa[j] != sa[p]) {
                        break;
                    }
                }
            }

            sa[pos] = len2 == 0 ? ~alph : ~++alph;

            ptr1 = ptr2;
        }

        if (n + n <= h) {
            // count and compact names
            for (int i = k - 1; i >= h; i--) {
                int ptr = sa[i];
                if (ptr < 0) {
                    sa[i] = 0;
                    ptr = ~ptr;
                    sa[--k] = ptr;
                    sa[n + ptr]++;
                }
            }
        } else {
            int idx = k;
            for (int i = k - 1; i >= h; i--) {
                final int ptr = sa[i];
                if (ptr < 0) {
                    sa[i] = 0;
                    sa[--idx] = ~ptr;
                }
            }
            for (int i = idx; i < k; i++) {
                sa[n + sa[i]]++;
            }
        }

        return ++alph == n ? 0 : alph;
    }

    private static int nameSubstrSA_saka_k(final int[] sa, int k, final int m, final int n) {

        final int h = k >> 1; // part by half        
        int alph = -1;

        // name sorted LMS substrings (bk - bucket start)
        for (int i = 0, ptr1 = 0, len1 = 0, pos1 = 0, bk = 0; i < n; i++) {
            final int ptr2 = sa[i];
            sa[i] = 0;
            final int pos2 = h + (ptr2 >> 1);
            int len2 = getNextLMS(sa, ptr2 + k, m - 1) - ptr2 - k + 1;
            if (len1 != len2) {
                len1 = len2;
            } else {
                for (int j = ptr1 + k, p = ptr2 + k; len2 > 0; j++, p++, len2--) {
                    if (sa[j] != sa[p]) {
                        break;
                    }
                }
            }

            if (len2 == 0) { // equal substrings
                sa[pos1] = ~bk;
                sa[pos2] = ~bk;
                if (i == bk + 1) {
                    sa[bk] = -2;
                    sa[i] = ~alph;
                } else {
                    sa[bk]--;
                }
            } else {
                sa[i] = ++alph;
                sa[pos2] = ~i;
                bk = i;
            }
            ptr1 = ptr2;
            pos1 = pos2;
        }
        
        return alph + 1;
    }
    
    private static void compactNames_sais(final int sa[], int k, final int n, final int l) {
        final int h = k >> 1; // part by half
        
        if (n + l <= h) {
            for (int i = k - 1; i >= h; i--) {
                int ptr = sa[i];
                if (ptr < 0) {
                    sa[i] = 0;
                    ptr = ~ptr;
                    final int p = sa[ptr] < 0 ? ~sa[ptr + 1] : sa[ptr];
                    sa[--k] = p;
                    sa[n + p]++;
                }
            }
        } else {
            int idx = k;
            for (int i = k - 1; i >= h; i--) {
                int ptr = sa[i];
                if (ptr < 0) {
                    sa[i] = 0;
                    ptr = ~ptr;
                    sa[--idx] = sa[ptr] < 0 ? ~sa[ptr + 1] : sa[ptr];
                }
            }

            for (int i = idx; i < k; i++) {
                sa[n + sa[i]]++;
            }
        }
    }
    
    /**
     * <p>
     * Compacts new string and updates SACA-K alphabet based on 'L' / 'S'.
     * <br/>
     * (Nong, G. 2013 5. NAMING SORTED LMS-SUBSTRINGS)
     * </p>
     * <pre>
     * ccccc - previously calculated ~(count) for every sssss character
     * 
     *        0            h            k   sa.length
     * sa[] = |ccccc       |       sssss|...|
     *                             |
     *                             return
     * </pre>
     * 
     * @param sa the array to work with.
     * @param k the last position to store the string.
     * 
     * @return the first string position
     */
    private static int compactNames_saka_k(final int sa[], int k) {

        final int h = k >> 1; // part by half
        
        boolean t = false;
        for (int i = k - 1, b0 = Integer.MIN_VALUE; i >= h; i--) {
            int b1 = sa[i];
            sa[i] = 0;
            if (b1 < 0) {
                b1 = ~b1;
                if (b1 > b0) {
                    t = false;
                } else if (b1 < b0) {
                    t = true;
                    if (sa[b1] < 0) {
                        b1 += ~sa[b1];
                    }
                } else if (t) {
                    if (sa[b1] < 0) {
                        b1 += ~sa[b1];
                    }
                }
                sa[--k] = b1;
                b0 = b1;
            }
        }
        return k;
    }
    
    /**
     * Finds the next LMS position in the original byte array
     * 
     * @param s the the original byte array
     * @param i the current 'S*' position
     * 
     * @return the next 'S*' position
     */
    private static int getNextLMS(final ByteBuffer s, int i) {
        final int m = s.limit() - 1;
        
        while(i < m && (s.get(i) & 0xFF) <= (s.get(++i) & 0xFF)); // find 'L'
        
        int pos;
        byte val;
        do {
            pos = i;
            val = s.get(i);
            while (i++ < m && val == s.get(i));
        } while (i <= m && (val & 0xFF) > (s.get(i) & 0xFF));

        return pos;
    }

    private static int getNextLMS(final int[] sa, int i, final int m) {
        
        while(i < m && sa[i] <= sa[++i]); // find 'L'
        
        int pos;
        int val;
        do {
            pos = i;
            val = sa[i];
            while (i++ < m && val == sa[i]);
        } while (i <= m && val > sa[i]);

        return pos;
    }
}
