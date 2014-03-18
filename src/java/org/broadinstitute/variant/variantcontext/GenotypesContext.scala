/*
* Copyright (c) 2012 The Broad Institute
*
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.variantcontext

import com.google.java.contract.Ensures
import com.google.java.contract.Requires
import org.broad.tribble.util.ParsingUtils
import scala.collection.mutable.ArrayBuffer
import scala.collection.Set


/**
 * Represents an ordered collection of Genotype objects
 */

/**
 * Ported to Scala by Wim Spee on 2/24/14.
 */
/**
 * Create a fully resolved GenotypeContext containing genotypes, sample lookup table,
 * and sorted sample names
 *
 * @param genotypes our genotypes in arbitrary
 * @param sampleNameToOffset map optimized for efficient lookup.  Each genotype in genotypes must have its
 *                           sample name in sampleNameToOffset, with a corresponding integer value that indicates the offset of that
 *                           genotype in the vector of genotypes
 * @param sampleNamesInOrder a list of sample names, one for each genotype in genotypes, sorted in alphabetical
 *                           order.
 */
@Requires(Array[String]("genotypes != null",
  "sampleNameToOffset != null",
  "sampleNamesInOrder != null",
  "genotypes.size() == sampleNameToOffset.size()",
  "genotypes.size() == sampleNamesInOrder.size()"))
class GenotypesContext protected(
                                  private val genotypes: ArrayBuffer[Genotype],
                                  private var sampleNameToOffset: scala.collection.mutable.Map[String, Integer], //a map optimized for efficient lookup.  Each genotype in genotypes must have its sample name in sampleNameToOffset, with a corresponding integer value that indicates the offset of that genotype in the vector of genotypes
                                  private var sampleNamesInOrder: Array[String] //sampleNamesInOrder a list of sample names, one for each genotype in genotypes, sorted in alphabetical order
                                  ) extends scala.collection.mutable.Buffer[Genotype] {


  // ---------------------------------------------------------------------------
  //
  // protected auxillary constructors -- you have to use static create methods to make these classes
  //
  // ---------------------------------------------------------------------------
  @Requires(Array[String]("genotypes != null"))
  protected def this(genotypes: ArrayBuffer[Genotype])      { this(genotypes, null, null)         }   //Create an GenotypeContext containing genotypes
  @Requires(Array[String]("n >= 0"))
  protected def this(n: Int)                                { this(new ArrayBuffer[Genotype](n)); }   //Create an empty GenotypeContext, with initial capacity for n elements
  protected def this()                                      { this(10); }                             //Create an empty GenotypeContext


  /**
   * An ArrayList of genotypes contained in this context
   *
   * WARNING: TO ENABLE THE LAZY VERSION OF THIS CLASS, NO METHODS SHOULD DIRECTLY
   * ACCESS THIS VARIABLE.  USE getGenotypes() INSTEAD.
   *
   */
  //val notToBeDirectlyAccessedGenotypes = genotypes;

  private var genotypeCounts : Array[Int] = null;      // Counts for each of the possible Genotype types in this context
  private var maxPloidy = -1;                          // Cached value of the maximum ploidy observed among all samples

  private var mNotYetComputed = true
  private var monomorphic : Boolean = false;


  // ---------------------------------------------------------------------------
  //
  // Monomorphic in samples methods
  //
  // ---------------------------------------------------------------------------


  /**
   * Genotype-specific functions -- are the genotypes monomorphic w.r.t. to the alleles segregating at this
   * site?  That is, is the number of alternate alleles among all fo the genotype == 0?
   *
   * @return true if it's monomorphic
   */
  def isMonomorphicInSamples(reference : Allele) : Boolean = {
    if ( mNotYetComputed ){
      monomorphic =  hasGenotypes && getCalledChrCount(reference) == getCalledChrCount();
      mNotYetComputed = false
    }
    monomorphic;
  }

  /**
   * Genotype-specific functions -- are the genotypes polymorphic w.r.t. to the alleles segregating at this
   * site?  That is, is the number of alternate alleles among all fo the genotype > 0?
   *
   * @return true if it's polymorphic
   */
  def isPolymorphicInSamples(reference : Allele ): Boolean = { ! isMonomorphicInSamples(reference);  }



  // ---------------------------------------------------------------------------
  //
  // Genotype counts methods
  //
  // ---------------------------------------------------------------------------

  /**
   * Count and store the amount of times that a genotype is NO_CALL, HOM_REF, HET, HOM_VAR or MIXED   *
   */
  private def calculateGenotypeCounts() {
    if ( genotypeCounts == null )
    {
      genotypeCounts = genotypes.foldLeft(new Array[Int](6)){(accum, element) => val idx = element.getType().ordinal() ; val oldValue = accum(idx); accum.update(idx, oldValue +1); accum  }
    }
  }
  /**
   * Genotype-specific functions -- how many no-calls are there in the genotypes?   *
   * @return number of no calls
   */
  def getNoCallCount() :Int  ={ calculateGenotypeCounts(); genotypeCounts(GenotypeType.NO_CALL.ordinal())  }
  def getHomRefCount() :Int  ={ calculateGenotypeCounts(); genotypeCounts(GenotypeType.HOM_REF.ordinal())  }
  def getHetCount() : Int    ={ calculateGenotypeCounts(); genotypeCounts(GenotypeType.HET.ordinal())      }
  def getHomVarCount() : Int ={ calculateGenotypeCounts(); genotypeCounts(GenotypeType.HOM_VAR.ordinal())  }
  def getMixedCount() : Int  ={ calculateGenotypeCounts(); genotypeCounts(GenotypeType.MIXED.ordinal())    }




  // ---------------------------------------------------------------------------
  //
  // Called chromosome  counts methods
  //
  // ---------------------------------------------------------------------------


  /**
   * Returns the number of chromosomes carrying any allele in the genotypes (i.e., excluding NO_CALLS)
   *
   * @return chromosome count
   */
  def getCalledChrCount() : Int = { getCalledChrCount(GenotypesContext.EMPTY_SET); }

  /**
   * Returns the number of chromosomes carrying any allele in the genotypes (i.e., excluding NO_CALLS)
   *
   * @param sampleIds IDs of samples to take into account. If empty then all samples are included.
   * @return chromosome count
   */
  def getCalledChrCount( sampleIds : scala.collection.Set[String]) : Int = {

    val genotypes = if(sampleIds.isEmpty){ this.genotypes}else{ getGenotypesByName(sampleIds);}
    genotypes.flatMap(_.getAlleles()).count(_.isCalled)
  }

  /**
   * Returns the number of chromosomes carrying allele A in the genotypes
   *
   * @param a allele
   * @return chromosome count
   */
  def getCalledChrCount( a : Allele)  : Int = {
    getCalledChrCount(a,GenotypesContext.EMPTY_SET);
  }

  /**
   * Returns the number of chromosomes carrying allele A in the genotypes
   *
   * @param a allele
   * @param sampleIds - IDs of samples to take into account. If empty then all samples are included.
   * @return chromosome count
   */
  def getCalledChrCount( a : Allele,  sampleIds : Set[String])  : Int = {

    val genotypes = if(sampleIds.isEmpty){ this.genotypes} else{ getGenotypesByName(sampleIds);}
    genotypes.map(_.countAllele(a)).sum
  }


  // ---------------------------------------------------------------------------
  //
  // MaxPloidy methods
  //
  // ---------------------------------------------------------------------------



  /**
   * Determine the max ploidy. 0 if genotypes is empty.
   */
  def determineMaxPloidy() =  if (genotypes.isEmpty){ maxPloidy = 0 }else{ maxPloidy = genotypes.map(_.getPloidy()).max }

  /**
   * What is the max ploidy among all samples?  Returns defaultPloidy if no genotypes are present
   *
   * @param defaultPloidy the default ploidy, if all samples are no-called
   * @return
   */
  @Ensures(Array[String]("result >= 0"))
  def getMaxPloidy( defaultPloidy :  Int) :Int  = if(maxPloidy == 0) defaultPloidy else maxPloidy


  // ---------------------------------------------------------------------------
  //
  // Mutability methods
  //
  // ---------------------------------------------------------------------------
  private var _immutable = false;     /** Are we allowing users to modify the list? */
  def immutable(): GenotypesContext = { _immutable = true;  this;   }
  def isMutable                     = { !_immutable;                }
  def checkImmutability()             { if (_immutable) { throw new IllegalAccessError("GenotypeMap is currently immutable, but a mutator method was invoked on it"); }  }

  // ---------------------------------------------------------------------------
  //
  // caches
  //
  // ---------------------------------------------------------------------------

  //invalidate chaches
  protected def invalidateCaches() { invalidateSampleNameMap(); invalidateSampleOrdering()}

  @Ensures(Array[String]("sampleNameToOffset == null"))
  protected def invalidateSampleNameMap()   {  sampleNameToOffset = null;  }

  @Ensures(Array[String]("sampleNamesInOrder == null"))
  protected def invalidateSampleOrdering()  {  sampleNamesInOrder = null;  }

  //create caches
  @Ensures(Array[String]( "sampleNamesInOrder != null" ))
  protected def ensureSampleOrdering() {
    if (sampleNamesInOrder == null) {   sampleNamesInOrder = genotypes.map(_.getSampleName()).toArray.sortWith(_ < _ ) }
  }

  @Ensures(Array[String]("sampleNameToOffset != null"))
  protected def ensureSampleNameMap() {
    if (sampleNameToOffset == null) {
      //create a new sampleNameToOffset HashMap with for every genotype a mapping of the sample name to the position in the arraybuffer.
      //position in the buffer is the same as the current size of the hashmap
      sampleNameToOffset = genotypes.foldLeft(scala.collection.mutable.Map[String,Integer]()){(accum, element) => accum + (element.getSampleName() -> accum.size)}
    }
  }



  def getSampleI(sampleName: String): Integer = { ensureSampleNameMap();  sampleNameToOffset.get(sampleName).get;   }

  /**
   * Iterate over the Genotypes in this context in the order specified by sampleNamesInOrder
   *
   * @param sampleNamesInOrder a Iterable of String, containing exactly one entry for each Genotype sample name in
   *                           this context
   * @return a Iterable over the genotypes in this context.
   */
  @Requires(Array[String]("sampleNamesInOrder != null"))
  def iterateInSampleNameOrder(sampleNamesInOrder: Iterable[String]): Iterable[Genotype] = {
    import scala.collection.JavaConversions._
    new Iterable[Genotype]()
    {
      override def iterator: Iterator[Genotype] = {  new InOrderIterator(sampleNamesInOrder.iterator);  }
    }
  }

  /**
   * Iterate over the Genotypes in this context in their sample name order (A, B, C)
   * regardless of the underlying order in the vector of genotypes
   * @return a Iterable over the genotypes in this context.
   */
  def iterateInSampleNameOrder(): Iterable[Genotype] = {
    import scala.collection.JavaConversions._
    iterateInSampleNameOrder(getSampleNamesOrderedByName())
  }

  private class InOrderIterator(private val sampleNamesInOrder: Iterator[String]) extends Iterator[Genotype] {

    override def hasNext: Boolean =    { sampleNamesInOrder.hasNext                 }
    override def next(): Genotype =    { apply(sampleNamesInOrder.next())           }
    def remove()              { throw new UnsupportedOperationException()  }
  }

  /**
   * @return The set of sample names for all genotypes in this context, in arbitrary order
   */
  @Ensures(Array[String]("result != null"))
  def getSampleNames(): scala.collection.Set[String] =
  {
    ensureSampleNameMap(); sampleNameToOffset.keySet
  }

  /**
   * @return The set of sample names for all genotypes in this context, in their natural ordering (A, B, C)
   */
  @Ensures(Array[String]("result != null"))
  def getSampleNamesOrderedByName(): Array[String] =
  {
    ensureSampleOrdering();  sampleNamesInOrder
  }

  @Requires(Array[String]("sample != null"))
  def containsSample(sample: String): Boolean =
  {
    ensureSampleNameMap(); sampleNameToOffset.contains(sample)
  }

  @Requires(Array[String]("samples != null"))
  def containsSamples(otherSamples: Iterable[String]): Boolean =
  {
    val thisSamples = getSampleNames()
    val countContained = otherSamples.count( name => thisSamples.contains(name))

    if(countContained == otherSamples.size){true}else{false}
  }


  @Requires(Array[String]("samples != null"))
  @Ensures(Array[String]("result != null"))
  def subsetToSamples(samples: Array[String]): GenotypesContext = {
    subsetToSamples(samples.toSet)
  }

  /**
   * Return a freshly allocated subcontext of this context containing only the samples
   * listed in samples.  Note that samples can contain names not in this context, they
   * will just be ignored.
   *
   * @param samples
   * @return
   */
  @Requires(Array[String]("samples != null"))
  @Ensures(Array[String]("result != null"))
  def subsetToSamples(samples: Set[String]): GenotypesContext = {
    if (samples.size == 0) {
      GenotypesContext.NO_GENOTYPES
    }
    else {
      val subSetGenotypes = genotypes.filter(g => samples.contains(g.getSampleName()))
      GenotypesContext.create(subSetGenotypes)
    }
  }

  @Requires(Array[String]("samples != null"))
  @Ensures(Array[String]("result != null"))
  def getGenotypesByName(samples: Seq[String]): Seq[Genotype] = {
    getGenotypesByName(samples.toSet)
  }


  @Requires(Array[String]("samples != null"))
  @Ensures(Array[String]("result != null"))
  def getGenotypesByName(samples: Set[String]): Seq[Genotype] = {
    if (samples.size == 0)  { GenotypesContext.NO_GENOTYPES_ARRAY                         }
    else                    { genotypes.filter(g => samples.contains(g.getSampleName()))  }
  }

  /**
   * @return true if the context has associated genotypes
   */
  def hasGenotypes: Boolean = { !genotypes.isEmpty }

  def getGenotypes : Array[Genotype] = {  genotypes.toArray     }

  // ---------------------------------------------------------------------------
  //
  // Lazy methods
  //
  // ---------------------------------------------------------------------------

  //    public boolean isLazyWithData() {
  //        return this instanceof LazyGenotypesContext &&
  //                ((LazyGenotypesContext)this).getUnparsedGenotypeData() != null;
  //    }

  // ---------------------------------------------------------------------------
  //
  //  simple mappings of buffer functions
  //
  // ---------------------------------------------------------------------------
  override def size                                           = { genotypes.size              }
  override def length                                         = { genotypes.length            }
  override def isEmpty                                        = { genotypes.isEmpty           }
  override def contains(o: Any)                               = { genotypes.contains(o)       }
  override def iterator: Iterator[Genotype]                   = { genotypes.iterator          }
  override def apply(i: Int): Genotype                        = { genotypes(i)                }
  override def slice(i: Int, i1: Int): ArrayBuffer[Genotype]  = { genotypes.slice(i, i1);     }
  override def +=:(g : Genotype)                              = {throw new UnsupportedOperationException(); this}


  /**
   * Gets sample associated with this sampleName, or null if none is found   *
   * @param sampleName
   * @return
   */
  def apply(sampleName: String): Genotype = {
    val offset = getSampleI(sampleName)
    if (offset == null) { null                        }
    else                { genotypes(offset) }
  }


  private def containsAny(genotypes: Array[_ <: Genotype]): Boolean = { val count = genotypes.count(contains(_)); if(count > 0){true}else{false} ; }

  //unsupported mappings
  @Requires(Array[String]("! contains(genotype)"))
  override def insert(i: Int, genotype: Genotype*)  = { throw new UnsupportedOperationException() }
  override def insertAll(i: Int, genotypes: collection.Traversable[Genotype])  { throw new UnsupportedOperationException()  }

  /**
   * Note that remove requires us to invalidate our sample -> index cache.
   *
   * Looping over multiple samples and calling the remove method below  is extremely inefficient,
   * as each call to remove invalidates the cache * and containsSample requires us to rebuild it, an O(n) operation.
   *
   * If you must remove many samples from the GC, use either removeAll or retainAll  * to avoid this O(n * m) operation.
   *
   * @param i
   * @return
   */
  override def remove(i: Int): Genotype =                   { checkImmutability(); invalidateCaches(); genotypes.remove(i)                    }
  override def -=(g: Genotype)  =                           { checkImmutability(); invalidateCaches(); genotypes -= g; this                   }
  override def --=(gArray: TraversableOnce[Genotype])  =    { checkImmutability(); invalidateCaches(); genotypes --= gArray; this             }
  override def clear()                                      { checkImmutability(); invalidateCaches(); genotypes.clear()                      }
  def retainAll(gT: Seq[Genotype]) =                        { checkImmutability(); invalidateCaches(); genotypes --= genotypes.diff(gT); this }



  /**
   * Adds a single genotype to this context.
   *
   * There are many constraints on this input, and important
   * impacts on the performance of other functions provided by this
   * context.
   *
   * First, the sample name of genotype must be unique within this
   * context.  However, this is not enforced in the code itself, through
   * you will invalid the contract on this context if you add duplicate
   * samples and are running with CoFoJa enabled.
   *
   * Second, adding genotype also updates the sample name -> index map,
   * so add() followed by containsSample and related function is an efficient
   * series of operations.
   *
   * Third, adding the genotype invalidates the sorted list of sample names, to
   * add() followed by any of the SampleNamesInOrder operations is inefficient, as
   * each SampleNamesInOrder must rebuild the sorted list of sample names at
   * an O(n log n) cost.
   *
   * @param genotype
   * @return
   */

  @Requires(Array[String]("genotype != null", "get(genotype.getSampleName()) == null"))
  override
  def +=(genotype: Genotype) = {
    checkImmutability()
    invalidateSampleOrdering()

    if (sampleNameToOffset != null) { sampleNameToOffset += (genotype.getSampleName()-> size)}    // update the name map by adding entries
    genotypes += genotype
    this
  }


  /**
   * Adds all of the genotypes to this context   *
   * See {@link #add(Genotype)} for important information about this functions
   * constraints and performance costs
   *
   * @param genotypes
   * @return
   */

  @Requires(Array[String]("! containsAny(genotypes)"))
  override
  def ++=(genotypes: TraversableOnce[_ <: Genotype]) = {
    checkImmutability()
    invalidateSampleOrdering()

    if (sampleNameToOffset != null) {
      var pos = size
      genotypes.foreach { g => sampleNameToOffset += (g.getSampleName() -> pos); pos += 1  }
    }
    this.genotypes ++= genotypes
    this
  }

  override
  def update(i: Int, genotype: Genotype) {
    checkImmutability()
    val prev = genotypes(i)
    genotypes.update(i, genotype)

    invalidateSampleOrdering()
    if (sampleNameToOffset != null) {
      // update the name map by removing the old entry and replacing it with the new one
      sampleNameToOffset - (prev.getSampleName())
      sampleNameToOffset += (genotype.getSampleName()-> i)
    }

  }

  /**
   * Replaces the genotype in this context -- note for efficiency
   * reasons we do not add the genotype if it's not present.  The
   * return value will be null indicating this happened.
   *
   * Note this operation is preserves the map cache Sample -> Offset but
   * invalidates the sorted list of samples.  Using replace within a loop
   * containing any of the SampleNameInOrder operation requires an O(n log n)
   * resorting after each replace operation.
   *
   * @param genotype a non null genotype to bind in this context
   * @return null if genotype was not added, otherwise returns the previous genotype
   */
  @Requires(Array[String]("genotype != null"))
  def replace(genotype: Genotype) {
    checkImmutability()
    val offset = getSampleI(genotype.getSampleName())
    if (offset != null) { update(offset, genotype) }
  }

  override def toString(): String = {
    val gsArray  : Array[String] = iterateInSampleNameOrder.map(_.toString()).toArray
    val string = "[" + ParsingUtils.join(",", gsArray) + "]"
    string
  }

  def isLazyWithData : Boolean =  { false }

  def getUnparsedGenotypeData : String = {""}


}

object GenotypesContext {

  import java.util._
  /**
   * static constant value for an empty GenotypesContext.  Useful since so many VariantContexts have no genotypes
   */
  val NO_GENOTYPES = new GenotypesContext( ArrayBuffer[Genotype](), scala.collection.mutable.Map[String, Integer](), Array[String]() ).immutable()
  val EMPTY_SET = scala.collection.immutable.HashSet[String]()
  val NO_GENOTYPES_ARRAY = Array[Genotype]()

  // ---------------------------------------------------------------------------
  //
  // public static factory methods
  //
  // ---------------------------------------------------------------------------

  /**
   * Basic creation routine
   * @return an empty, mutable GenotypeContext
   */
  @Ensures(Array[String]("result != null"))
  def create(): GenotypesContext =               {   new GenotypesContext();  }

  /**
   * Basic creation routine
   * @return an empty, mutable GenotypeContext with initial capacity for nGenotypes
   */
  @Requires(Array[String]("nGenotypes >= 0"))
  @Ensures(Array[String]("result != null"))
  def create(nGenotypes: Int): GenotypesContext = {  new GenotypesContext(nGenotypes);  }

  /**
   * Create a fully resolved GenotypeContext containing genotypes, sample lookup table,
   * and sorted sample names
   *
   * @param genotypes our genotypes in arbitrary
   * @param sampleNameToOffset map optimized for efficient lookup.  Each genotype in genotypes must have its
   *                           sample name in sampleNameToOffset, with a corresponding integer value that indicates the offset of that
   *                           genotype in the vector of genotypes
   * @param sampleNamesInOrder a list of sample names, one for each genotype in genotypes, sorted in alphabetical
   *                           order.
   * @return an mutable GenotypeContext containing genotypes with already present lookup data
   */
  @Requires(Array[String]("genotypes != null", "sampleNameToOffset != null", "sampleNamesInOrder != null"))
  @Ensures(Array[String]("result != null"))
  def create(genotypes: ArrayBuffer[Genotype], sampleNameToOffset: scala.collection.mutable.Map[String, Integer], sampleNamesInOrder: Array[String]): GenotypesContext = {
    new GenotypesContext(genotypes, sampleNameToOffset, sampleNamesInOrder)
  }

  def create(genotypes: Array[Genotype], sampleNameToOffset: scala.collection.mutable.Map[String, Integer], sampleNamesInOrder: Array[String]): GenotypesContext = {
    val buffer = scala.collection.mutable.ArrayBuffer[Genotype]()
    genotypes.foreach(buffer += _)

    create(buffer, sampleNameToOffset, sampleNamesInOrder)
  }

  /**
   * Create a fully resolved GenotypeContext containing genotypes
   *
   * @param genotypes our genotypes in arbitrary
   * @return an mutable GenotypeContext containing genotypes
   */
  @Requires(Array[String]("genotypes != null"))
  @Ensures(Array[String]("result != null"))
  def create(genotypes: ArrayBuffer[Genotype]): GenotypesContext = {
    if (genotypes == null) { NO_GENOTYPES                   }
    else                   { new GenotypesContext(genotypes)}
  }

  /**
   * Create a fully resolved GenotypeContext containing genotypes
   *
   * @param genotypes our genotypes in arbitrary
   * @return an mutable GenotypeContext containing genotypes
   */
  @Requires(Array[String]( "genotypes != null"))
  @Ensures(Array[String]( "result != null"))
  def create(genotypes: Genotype*): GenotypesContext =
  { val arrayBuffer = ArrayBuffer[Genotype]()
    genotypes.foreach( arrayBuffer += _)
    create(arrayBuffer)  }

  /**
   * Create a freshly allocated GenotypeContext containing the genotypes in toCopy
   *
   * @param toCopy the GenotypesContext to copy
   * @return an mutable GenotypeContext containing genotypes
   */
  @Requires(Array[String]("toCopy != null"))
  @Ensures(Array[String]("result != null"))
  def copy(toCopy: GenotypesContext): GenotypesContext =
  {
    val copiedArrayBuffer = ArrayBuffer[Genotype]()
    toCopy.genotypes.copyToBuffer(copiedArrayBuffer)
    create(copiedArrayBuffer)
  }

  /**
   * Create a GenotypesContext containing the genotypes in iteration order contained
   * in toCopy
   *
   * @param toCopy the collection of genotypes
   * @return an mutable GenotypeContext containing genotypes
   */
  @Ensures(Array[String]("result != null"))
  def copy(toCopy: Iterable[Genotype]): GenotypesContext = {
    if  (toCopy == null)
    { NO_GENOTYPES }
    else
    {
      val copiedArrayBuffer = ArrayBuffer[Genotype]()
      toCopy.copyToBuffer(copiedArrayBuffer)
      create( copiedArrayBuffer)
    }
  }

}



