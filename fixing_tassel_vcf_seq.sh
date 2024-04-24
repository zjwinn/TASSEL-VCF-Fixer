#!/bin/bash

# Function to display usage instructions
usage() {
    echo
    echo "###############################################"
    echo "#                                             #"
    echo "#               Help Information              #"
    echo "#                                             #"
    echo "###############################################"
    echo
    echo "Usage:"
    echo -e "\t$0 [OPTIONS] ARGUMENT"
    echo
    echo "Description:"
    echo -e "\tThis script will take a known reference genome (.fa) and a molecular marker file (.vcf.gz)"
    echo -e "\tand checks for agreement between the marker file and the reference sequences. This function was made in"
    echo -e "\treponse to the TASSEL v2.0 Standalone GBS pipeline only putting VCFs out in major/minor allele format."
    echo -e "\tThis script ensures that the VCF shares the same ref/alt information as the reference genome."
    echo
    echo "Options:"
    echo -e "\t-h, --help               Display this help message and exit."
    echo
    echo "Arguments:"
    echo -e "\t-v, --vcf,               The VCF input (.vcf.gz file) file path. This must be the full length real path!"
    echo -e "\t-r, --ref,               The reference sequence (.fa file) file path that you are checking against the VCF. This must be the full length real path!"
    echo -e "\t-n, --name               The name of the new vcf.gz output."
    echo     
    echo "Examples:"
    echo -e "\tbash $0 -v 'example.vcf.gz' -r 'example.fa' -n 'fixed_example.vcf.gz'"
    echo
    exit 1
}

# Initialize variables
vcf=""
refseq=""
cores=""
name=""

# Parse command-line options
while getopts "v:r:n:h" opt; do
    case ${opt} in
        v | --vcf )
            vcf=$OPTARG
            ;;
        r | --ref )
            refseq=$OPTARG
            ;;    
        n | --name )    
            name=$OPTARG
            ;;
        h | --help )
            usage
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done
shift $((OPTIND -1))

# Check if vcf and refseq are provided
if [ -z "$vcf" ] || [ -z "$refseq" ]; then
    echo
    echo "*** Error: Both VCF and reference sequence files are required! ***"
    usage
fi

# Check if the files exist
if [ ! -f "$vcf" ]; then
    echo
    echo "*** Error: VCF file '$vcf' not found! ***"
    usage
fi

if [ ! -f "$refseq" ]; then
    echo
    echo "*** Error: Reference sequence '$refseq' not found! ***"
    usage
fi

# Echo header
echo 
echo "###############################################"
echo "#                                             #"
echo "#          TASSEL VCF Fixer v1.0              #"
echo "#                                             #"
echo "###############################################"
echo
echo "Written by: Zachary J. Winn PhD"
echo "Contact information:"
echo -e "\tGovernment Email: zachary.winn@usda.gov"
echo -e "\tPersonal Email: zwinn@outlook.com"
echo
echo "###############################################"
echo "# WARNING: This program is not under warranty #"
echo "#          Use at your own discretion!        #"
echo "###############################################"
echo

# Function to check if a string is a number
is_number() {
    # Regex pattern to match a number
    number_pattern="^[0-9]+([.][0-9]+)?$"
    
    # Check if the input matches the number pattern
    if [[ $1 =~ $number_pattern ]]; then
        echo 0
    else
        echo 1
    fi
}

# Main script
echo "### Provided path of VCF.GZ file: $vcf"
echo
echo "### Provided path of reference sequence FA file: $refseq"
echo 

# Print
echo "### Copying, unzipping, and editing vcf file..."
echo

# Make Temp Dir
if [ -d "fixing_tassel_vcf_temp" ]; then
    echo "*** Error: The directory "fixing_tassel_vcf_temp" already exist! Delete that directory and try again! ***"
    exit 1
fi

# Make a temp directory
mkdir fixing_tassel_vcf_temp

# Write out temp file which has the chromsomes listed as chr1A, chr1B, etc
cd ${PWD}/fixing_tassel_vcf_temp
cp $vcf ${PWD}/temp1.vcf.gz
gunzip temp1.vcf.gz

# Change id in temp1.vcf.gz if needed
if ! grep -q "##contig=<ID=Chr" temp1.vcf; then

    echo "*** Warning: string '##contig=<ID=Chr' not found in header, attempting to fix..."
    echo ""
    sed -i 's/##contig=<ID=/##contig=<ID=Chr/g' temp1.vcf

fi

if grep -q "UNKNOWN" temp1.vcf; then

    echo "*** Warning: string 'UNKNOWN' found in header, attempting to fix..."
    echo ""
    sed -i 's/UNKNOWN/Unknown/g' temp1.vcf

fi

# Add Chr to each chromosome
awk -v OFS='\t' '{ if ($1 !~ /^#/) $1 = "Chr" $1; print }' temp1.vcf > temp_file && mv temp_file temp1.vcf

# Make into GZ
bgzip temp1.vcf

# Print
echo "### Pulling positions from sequences for the provided VCF..."
echo

# Make a correpsonding list of the correct ref and alt in the reference
bedtools getfasta -fi $refseq -bed temp1.vcf.gz > correct_ref_alt.txt

# Print
echo "### Looping through and fixing ref/alt marker info..."
echo 

# Separate vcf header
zcat temp1.vcf.gz | grep --text '^#' > header.txt

# Separate vcf body
zcat temp1.vcf.gz | grep -v '^#' > body.txt

# Get marker info                                                                                                                                                                                                                                                                                                                              
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' "temp1.vcf.gz" > marker_info.txt

# Find length
length=$(zcat temp1.vcf.gz | grep -v '^#' | wc -l)

# Make an empty body
touch fixed_body.vcf

# Loop and fix
for i in $(seq 1 $length); do

    # Pull
    CHROM=$(awk -F'\t' -v row=$i -v col=1 'NR==row {print $col}' "marker_info.txt")
    POS=$(awk -F'\t' -v row=$i -v col=2 'NR==row {print $col}' "marker_info.txt")
    BACKPOS=$((POS-1))
    ID=$(awk -F'\t' -v row=$i -v col=3 'NR==row {print $col}' "marker_info.txt")
    REF=$(awk -F'\t' -v row=$i -v col=4 'NR==row {print $col}' "marker_info.txt")
    ALT=$(awk -F'\t' -v row=$i -v col=5 'NR==row {print $col}' "marker_info.txt")
    BEDREF=$(echo ">${CHROM}:${BACKPOS}-${POS}")
    BEDLOC=$(grep -n "$BEDREF" "correct_ref_alt.txt" | cut -d: -f1)
    TRUEREF=$((BEDLOC+1))
    TRUEREF=$(awk "NR == $TRUEREF" "correct_ref_alt.txt")

    # Check if the GBS Ref and True Ref agree
    if [ "$REF" = "$TRUEREF" ]; then

        # Pull that line and put it in the new vcf
        NEWLINE=$(awk -F'\t' "NR == $i" "body.txt")
        echo "$NEWLINE" >> fixed_body.vcf

    elif [ "$REF" != "$TRUEREF" ] && [ "$ALT" != "$TRUEREF" ]; then

        # Pull that line, fix the ref alt, and then put in the new vcf
        echo "********************************************************************"
        echo "Error:"
        echo "${ID} (Line = ${i}) is not in correct Ref/Alt!"
        echo "The Reference provided states ${TRUEREF} is the correct allele..."
        echo "However, the GBS provided states that Ref=${REF} and Alt=${ALT}..."
        echo "Something is wrong with this marker! Omitting from the final VCF..."
        echo "********************************************************************"

    else

        #do something else
        echo "${ID} (Line = ${i}) is not in correct Ref/Alt! Fixing!"
        NEWLINE=$(awk -F'\t' "NR == $i" "body.txt")
        NEWLINE=($NEWLINE)
        NEWLINE[3]=$TRUEREF
        NEWLINE[4]=$REF
        NEWLINE_SPACED=$(echo "${NEWLINE[@]}" | tr ' ' '\t')
        echo "${NEWLINE_SPACED[@]}" >> fixed_body.vcf

    fi

done

# Bind header and body
cat "header.txt" >> "corrected_genotyping_file.vcf"
cat "fixed_body.vcf" >> "corrected_genotyping_file.vcf"

# Change id in temp1.vcf.gz
sed -i 's/##contig=<ID=Chr/##contig=<ID=/g' corrected_genotyping_file.vcf
sed -i 's/Unknown/UNKNOWN/g' corrected_genotyping_file.vcf
sed -i 's/Chr//g' corrected_genotyping_file.vcf
sed -i 's/Reference allele is not known. The major allele was used as reference allele/Alleles reported as reference-alternative/g' corrected_genotyping_file.vcf

# Zip it up using bgzip
bgzip corrected_genotyping_file.vcf

# Move to top directory
mv corrected_genotyping_file.vcf.gz ../

# Change directory back to top
cd ..

# Change name
mv corrected_genotyping_file.vcf.gz $name

# Index the vcf.gz
bcftools index -c $name

# Print 
echo "### Done!"

# Remove
rm -rf fixing_tassel_vcf_temp
