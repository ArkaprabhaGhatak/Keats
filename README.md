Keats
=====

Keats is a Scala port of the Picard library of the Broad Institute of Harvard and MIT.

The goal is to use the advanced language features of Scala to make the library easier to understand / use and more scalabe.
See Coursera Functional Programming Principles in Scala for the advanced Scala features. 
https://www.coursera.org/course/progfun
And of course Akka http://akka.io/

I already ported the core genomic data representation classes : VariantContext, Allele(Context) and Genotype(Context).
Porting can be done class by class and existing unit tests arereused. This is possible because Java and Scala are compatible.
I am now working on making the vcf/bcf codec, reader and writer threadsafe.

Keats is a work in progress and not ready for use. The code is MIT license. 
Feel free to fork (parts of) the project or help with porting the code to Scala.
