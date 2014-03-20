Keats
=====

Keats is a Scala port of the Picard library of the Broad Institute of Harvard and MIT. http://picard.sourceforge.net/

The goal is to use the advanced language features of Scala to make the library easier to understand / use and more scalabe for use in Genomics research. A side goal is to increase my knowledge of Scala and advanced representation of Genomic variation in data structures. 

See Coursera Functional Programming Principles in Scala for the advanced Scala features. 
https://www.coursera.org/course/progfun

And of course Akka http://akka.io/

I already ported the core genomic data representation classes : VariantContext, Allele(Context) and Genotype(Context).
Porting can be done class by class and existing unit tests are reused. This is possible because Java and Scala are compatible.
I am now working on making the vcf/bcf codec, reader and writer threadsafe.

Keats is a work in progress and not ready for use. The code is MIT license. 
Feel free to fork (parts of) the project or help with porting the code to Scala.
