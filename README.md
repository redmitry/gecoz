### **GE**nome **C**ompression **Z** - **GECOZ** is a set of various NGS formats and algorithms.

At the moment it contains an efficient FM-Index implementation.  

To build the tool:
```
gecoz/java/> mvn install
```
This will build a **gecoz/java/nova-gecoz/target/gecotools.jar** file.
___
```
gecoz/java/nova-gecoz/target/> java -jar gecotools.jar -h
```
Will get the command line help:
```
gecotools -i file [optional params]

parameters:

-h (--help)           - this help message
-i (--input)          - either *.fa or *.gcz
-o [header][from][to] - depends on the input parameters
                        (*.fa -> *.gcz, *.gcz -> *.fa, *gcz -> *.seq)
-c [header] 'string'  - count string occurrences in the *.gcz file
-s [header] 'string'  - search string in the *.gcz file
-t                    - use n threads
-v [level]            - verbose (default = WARNING)

examples:

>java -jar gecotools.jar -t 8 -i hg38.fa -o hg38.gcz
>java -jar gecotools.jar -i hg38.gcz -o hg38.fasta
>java -jar gecotools.jar -i hg38.gcz -o chr15.seq chr15
>java -jar gecotools.jar -i hg38.gcz -c ATTAACCCATGAAAA
>java -jar gecotools.jar -i hg38.gcz -s ATTAACCCATGAAAA
>java -jar gecotools.jar -i hg38.gcz -s chr11 ATTAACCCATGAAAA
```
