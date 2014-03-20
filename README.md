Keats
=====

Keats is a Scala port of the Picard library of the Broad Institute of Harvard and MIT. 
Picard is one of the most advanced and widely used libraries for representation and processing of Genomic data from BAM, VCF and BCF files. 
http://picard.sourceforge.net/

The goal of the port is to use the advanced language features of Scala to make the library easier to understand / use and more scalabe for use in Genomics research. A personal side goal is to increase my knowledge of Scala and advanced data structure representation of Genomic variation.

See Coursera Functional Programming Principles in Scala for the advanced Scala features. 
https://www.coursera.org/course/progfun

And of course Akka http://akka.io/

I already ported the core genomic data representation classes : VariantContext, Allele(Context) and Genotype(Context).
Porting can be done class by class and existing unit tests are reused. This is possible because Java and Scala are compatible.
I am now working on making the vcf/bcf codec, reader and writer threadsafe.

Keats is not yet end user ready software, but aimed at developers interested in Scala and the functionality of Picard.  The code is MIT license. 
Feel free to fork (parts of) the project or help with porting the code to Scala.
