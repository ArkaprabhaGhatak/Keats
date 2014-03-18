package org.broadinstitute.variant.variantcontext

/**
 * Created by Wim Spee on 2/27/14.
 */
class AlleleContext(val _alleles : Array[Allele]) extends Seq[Allele]
{


  /** A set of the alleles segregating in this context */
  // we need to make this a LinkedHashSet in case the user prefers a given ordering of alleles
  if ( _alleles == null ) { throw new IllegalArgumentException("Alleles cannot be null"); }
  protected val alleles : Array[Allele] = AlleleContext.makeAlleles(_alleles);

  // a fast cached access point to the ref / alt alleles for biallelic case
  private val REF : Allele = _alleles.filter(_.isReference)(0)

  // set to the alt allele when biallelic, otherwise == null
  private val ALT : Allele = if(_alleles.size ==2){_alleles.filter(!_.isReference)(0)}else{null}


  /** The type (cached for performance reasons) of this context */
  protected var variantType  : VariantType = null;



  // ---------------------------------------------------------------------------------------------------------
  //
  // Working with alleles
  //
  // ---------------------------------------------------------------------------------------------------------

  /**
   * @return the reference allele for this context
   */
  def getReference() = {
    if ( REF == null ){ throw new IllegalStateException("BUG: no reference allele found at " + this);}
    REF;
  }


  /**
   * @return true if the context is strictly bi-allelic
   */
  def isBiallelic() = { size == 2; }


  /**
   * @return The allele sharing the same bases as this String.  A convenience method; better to use byte[]
   */
  def  getAllele( allele : String) : Allele= { getAllele(allele.getBytes());  }

  /**
   * @return The allele sharing the same bases as this byte[], or null if no such allele is present.
   */
  def getAllele( allele : Array[Byte]) : Allele = { Allele.getMatchingAllele(getAlleles(), allele); }

  /**
   * @return True if this context contains Allele allele, or false otherwise
   */
  def hasAllele( allele : Allele) : Boolean =                                       { hasAllele(allele, false, true);           }
  def hasAllele( allele : Allele, ignoreRefState : Boolean) : Boolean =             { hasAllele(allele, ignoreRefState, true);  }
  def hasAlternateAllele( allele : Allele) : Boolean =                              { hasAllele(allele, false, false);          }
  def hasAlternateAllele( allele :Allele , ignoreRefState : Boolean) : Boolean  =   { hasAllele(allele, ignoreRefState, false); }

  private def hasAllele( allele : Allele, ignoreRefState : Boolean, considerRefAllele : Boolean) : Boolean =
  {
    if ( (considerRefAllele && allele == REF) || allele == ALT ) { true} // optimization for cached cases
    else
    {
      val allelesToConsider = if(considerRefAllele){ getAlleles()} else{ getAlternateAlleles();}
      val countEquals = allelesToConsider.count(_.equals(allele, ignoreRefState))
      if(countEquals > 0){true}else{false}
    }
  }

  def getNAlleles() = { size }

  /**
   * Gets the alleles.  This method should return all of the alleles present at the location,
   * including the reference allele.  There are no constraints imposed on the ordering of alleles
   * in the set. If the reference is not an allele in this context it will not be included.
   *
   * @return the set of alleles
   */
  def getAlleles() = {  alleles; }

  /**
   * Gets the alternate alleles.  This method should return all the alleles present at the location,
   * NOT including the reference allele.  There are no constraints imposed on the ordering of alleles
   * in the set.
   *
   * @return the set of alternate alleles
   */
  def getAlternateAlleles() =  { alleles.drop(1); }






  /**
   * @param i -- the ith allele (from 0 to n - 2 for a context with n alleles including a reference allele)
   * @return the ith non-reference allele in this context
   * @throws IllegalArgumentException if i is invalid
   */
  def getAlternateAllele( i: Int) ={ alleles(i+1);  }

  /**
   * @param  other  VariantContext whose alleles to compare against
   * @return true if this VariantContext has the same alleles (both ref and alts) as other,
   *         regardless of ordering. Otherwise returns false.
   */
  def hasSameAllelesAs (  other : AlleleContext ) : Boolean =  {
    if(hasSameAlternateAllelesAs(other) && other.getReference().equals(getReference(), false)) {true}
    else{false}
  }

  /**
   * @param  other  VariantContext whose alternate alleles to compare against
   * @return true if this VariantContext has the same alternate alleles as other,
   *         regardless of ordering. Otherwise returns false.
   */
  def hasSameAlternateAllelesAs (  other : AlleleContext) = {
    val thisAlternateAlleles = getAlternateAlleles();
    val otherAlternateAlleles = other.getAlternateAlleles();

    if ( thisAlternateAlleles.size != otherAlternateAlleles.size ) { false}
    else
    {
      val countAlleleContainedInOtherAlleles = thisAlternateAlleles.count(otherAlternateAlleles.contains(_))
      if(countAlleleContainedInOtherAlleles == thisAlternateAlleles.size){true}else{false}
    }
  }

  /**
   * Lookup the index of allele in this variant context
   *
   * @param allele the allele whose index we want to get
   * @return the index of the allele into getAlleles(), or -1 if it cannot be found
   */
  def getAlleleIndex(  allele : Allele) : Int  = { alleles.indexOf(allele) }

  /**
   * Return the allele index #getAlleleIndex for each allele in alleles
   *
   * @param alleles the alleles we want to look up
   * @return a list of indices for each allele, in order
   */
  def getAlleleIndices( alleles : Array[Allele]) : Array[Int] = { alleles.map(getAlleleIndex(_)) }

  def getGLIndecesOfAlternateAllele( targetAllele : Allele)  : Array[Int] = {
    val index = getAlleleIndex(targetAllele);
    if ( index == -1 ) throw new IllegalArgumentException("Allele " + targetAllele + " not in this VariantContex " + this);
    GenotypeLikelihoods.getPLIndecesOfAlleles(0, index);
  }

  def hasSymbolicAlleles  : Boolean = {  AlleleContext.hasSymbolicAllelesInArray(getAlleles());  }



  // ---------------------------------------------------------------------------
  //
  //  simple mappings of Seq functions
  //
  // ---------------------------------------------------------------------------

  override def size                                           = { alleles.size              }
  override def length                                         = { alleles.length            }
  override def iterator: Iterator[Allele]                     = { alleles.iterator          }
  override def apply(i: Int): Allele                          = { alleles(i)                }


  // ---------------------------------------------------------------------------
  //
  //  VariantType methods
  //
  // ---------------------------------------------------------------------------

  /**
   * Determines (if necessary) and returns the type of this variation by examining the alleles it contains.
   *
   * @return the type of this VariantContext
   **/
  def getType() : VariantType = { if ( variantType == null ){ determineType()}; variantType    }
  def isSNP()                       = { getType() == VariantType.SNP                                 }
  def isVariant()                   = { getType() != VariantType.NO_VARIATION                        }
  def isPointEvent()                = { isSNP() || !isVariant()                                      }
  def isIndel()                     = { getType() == VariantType.INDEL                               }
  def isSimpleInsertion()           = { isSimpleIndel() && getReference().length() == 1              } // can't just call !isSimpleDeletion() because of complex indels
  def isSimpleDeletion()            = { isSimpleIndel() && getAlternateAllele(0).length() == 1       } // can't just call !isSimpleInsertion() because of complex indels
  def isComplexIndel()              = { isIndel() && !isSimpleDeletion() && !isSimpleInsertion()     }
  def isSymbolic()                  = { getType() == VariantType.SYMBOLIC                            }
  def isSymbolicOrSV()              = { isSymbolic() || isStructuralIndel()                          }
  def isMNP()                       = { getType() == VariantType.MNP;                                }
  def isMixed()                     = { getType() == VariantType.MIXED;                              }

  /**
   * @return true if the alleles indicate a simple indel, false otherwise.
   */
  def isSimpleIndel() = {
    getType() == VariantType.INDEL  &&                 // allelic lengths differ
      isBiallelic()   &&                      // exactly 2 alleles
      getReference().length() > 0   &&        // ref is not null or symbolic
      getAlternateAllele(0).length() > 0  &&  // alt is not null or symbolic
      getReference().getBases()(0) == getAlternateAllele(0).getBases()(0)  &&  // leading bases match for both alleles
      (getReference().length() == 1 || getAlternateAllele(0).length() == 1);
  }

  def isStructuralIndel() : Boolean =  {
    var isStructuralIndelBoolean = false
    if ( getType() == VariantType.INDEL ) {
      val sizes = getIndelLengths();
      if ( sizes != null ) {
        val countTooLargeIndel = sizes.count(_ > AlleleContext.MAX_ALLELE_SIZE_FOR_NON_SV)
        if (countTooLargeIndel > 0){isStructuralIndelBoolean = true }
      }
    }
    isStructuralIndelBoolean;
  }

  /**
   * Gets the sizes of the alternate alleles if they are insertion/deletion events, and returns a list of their sizes
   *
   * @return a list of indel lengths ( null if not of type indel or mixed )
   */
  def getIndelLengths(): Array[Int] = {
    if (getType() != VariantType.INDEL && getType() != VariantType.MIXED) {
      null;
    }
    else {
      val refAlleleLenght = getReference().length()
      getAlternateAlleles().map(_.length() - refAlleleLenght)
    }
  }

  private def determineType() {
    if ( variantType == null ) {
      size match
      {
        case 0 => throw new IllegalStateException("Unexpected error: requested type of AlleleContext with no alleles!" + this);
        case 1 => variantType =  VariantType.NO_VARIATION;
        case _ => determinePolymorphicType();
      }
    }
  }

  private def determinePolymorphicType() {
    variantType = null;

    // do a pairwise comparison of all alleles against the reference allele
    var counter = 0
    val limit = size
    var mixedFound = false

    while(counter < limit && ! mixedFound)
    {
      val allele = alleles(counter)
      if ( allele == REF ){}
      else
      {
        // find the type of this allele relative to the reference
        val biallelicType = AlleleContext.typeOfBiallelicVariant(REF, allele);

        // for the first alternate allele, set the type to be that one
        if ( variantType == null ) {
          variantType = biallelicType;
        }
        // if the type of this allele is different from that of a previous one, assign it the MIXED type and quit
        else if ( biallelicType != variantType )
        {
          variantType = VariantType.MIXED;
          mixedFound = true
        }
      }
      counter +=1
    }
  }


  def validateAlleles() {

    val countRef = alleles.count(_.isReference)
    if(countRef == 0) { throw new IllegalArgumentException("No reference allele found in VariantContext");}
    if(countRef > 1)  { throw new IllegalArgumentException("BUG: Received two reference tagged alleles in VariantContext " + alleles + " this=" + this)}

    val countNoCall = alleles.count(_.isNoCall)
    if(countNoCall > 0 ){ throw new IllegalArgumentException("BUG: Cannot add a no call allele to a variant context " + alleles + " this=" + this);}
  }

}


object AlleleContext
{

  private val MAX_ALLELE_SIZE_FOR_NON_SV = 150;

  def hasSymbolicAllelesInArray( alleles : Array[Allele] ) : Boolean = { alleles.count(_.isSymbolic) > 0 }

  private def typeOfBiallelicVariant( ref : Allele, allele : Allele) : VariantType = {
    if ( ref.isSymbolic ){ throw new IllegalStateException("Unexpected error: encountered a record with a symbolic reference allele");}

    if ( allele.isSymbolic ) {VariantType.SYMBOLIC;}
    else
    {
      if ( ref.length() == allele.length() )
      {
        if ( allele.length() == 1 ){ VariantType.SNP;}
        else{ VariantType.MNP;}
      }
      else
      {
        // Important note: previously we were checking that one allele is the prefix of the other.  However, that's not an
        // appropriate check as can be seen from the following example:
        // REF = CTTA and ALT = C,CT,CA
        // This should be assigned the INDEL type but was being marked as a MIXED type because of the prefix check.
        // In truth, it should be absolutely impossible to return a MIXED type from this method because it simply
        // performs a pairwise comparison of a single alternate allele against the reference allele (whereas the MIXED type
        // is reserved for cases of multiple alternate alleles of different types).  Therefore, if we've reached this point
        // in the code (so we're not a SNP, MNP, or symbolic allele), we absolutely must be an INDEL.

        VariantType.INDEL;

        // old incorrect logic:
        // if (oneIsPrefixOfOther(ref, allele))
        //     return Type.INDEL;
        // else
        //     return Type.MIXED;
      }
    }
  }



  // protected basic manipulation routines
  private def makeAlleles( alleles : Array[Allele]) : Array[Allele] =
  {

    if ( alleles.isEmpty ){  throw new IllegalArgumentException("Cannot create a AlleleContext with an empty allele list");}

    val equalCounts = alleles.foldLeft(0){(accum, element) => accum + alleles.count(_.equals(element, true)) }

    if(equalCounts != alleles.size) {throw new IllegalArgumentException("Duplicate allele added to VariantContext: " + alleles)}

    val referenceCount = alleles.count(_.isReference)

    if(referenceCount == 0){ throw new IllegalArgumentException("Alleles for a AlleleContext must contain at least one reference allele: " + alleles) }

    if(referenceCount > 1){ throw new IllegalArgumentException("Alleles for a VariantContext must contain at most one reference allele: " + alleles)  }

    //put the reference allele in the first position of the array
    if(!alleles(0).isReference)
    {
       alleles.find(_.isReference).get +:  alleles.filterNot(_.isReference)
    }
    else
    {
      alleles
    }

  }


}






