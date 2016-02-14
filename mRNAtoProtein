import java.util.Scanner;

public class mRNAtoProteins {

	/**
	 * @description:
	 *	Get from the user input of mRNA. Should contains only rna codons as A, U, G, C and written fron 5' to 3' edges.
	 *	Prints the suitable protein (by amino acid).
	 * @throws:
	 *	IllegalArgumentException if the mRNA does not contains start codon (AUG) or end codon (UAA/UAG/UGA) or wrong bases
	 */
	public static void main(String[] args) {
		String mRNA;
		String protein = "";
		System.out.println("Enter mRNA from 5' to 3': ");
		Scanner sc = new Scanner(System.in);
		mRNA = sc.next().toUpperCase();

		for (int i=start(mRNA);i<end(mRNA);i=i+3){
			String codon = mRNA.substring(i, i+3);

			if (codon.equals("AUG"))
				protein=protein+"Met-";

			else if (codon.equals("UUU")||codon.equals("UUC"))
				protein=protein+"Phe-";

			else if ((codon.charAt(0)=='C' && codon.charAt(1)=='U')||
					(codon.equals("UUA")||codon.equals("UUG")))
				protein=protein+"Leu-";

			else if ((codon.charAt(0)=='U' && codon.charAt(1)=='C')||
					(codon.equals("AGU")||codon.equals("AGC")))
				protein=protein+"Ser-";

			else if (codon.equals("UAU")||codon.equals("UAC"))
				protein=protein+"Tyr-";

			else if (codon.equals("UGU")||codon.equals("UGC"))
				protein=protein+"Cyc-";

			else if (codon.equals("UGG"))
				protein=protein+"Trp-";

			else if (codon.charAt(0)=='C' && codon.charAt(1)=='C')
				protein=protein+"Pro-";

			else if (codon.equals("CAU")||codon.equals("CAC"))
				protein=protein+"His-";

			else if (codon.equals("CAA")||codon.equals("CAG"))
				protein=protein+"Gln-";

			else if ((codon.charAt(0)=='C' && codon.charAt(1)=='G')||
					(codon.equals("AGA")||codon.equals("AGG")))
				protein=protein+"Arg-";

			else if (codon.equals("AUU")||codon.equals("AUC")||codon.equals("AUA"))
				protein=protein+"Ile-";

			else if ((codon.charAt(0)=='A' && codon.charAt(1)=='C'))
				protein=protein+"Thr-";

			else if (codon.equals("AAU")||codon.equals("AAC"))
				protein=protein+"Asn-";

			else if (codon.equals("AAA")||codon.equals("AAG"))
				protein=protein+"Lys-";

			else if ((codon.charAt(0)=='G' && codon.charAt(1)=='U'))
				protein=protein+"Val-";

			else if ((codon.charAt(0)=='G' && codon.charAt(1)=='C'))
				protein=protein+"Ala-";

			else if (codon.equals("GAU")||codon.equals("GAC"))
				protein=protein+"Asp-";

			else if (codon.equals("GAA")||codon.equals("GAG"))
				protein=protein+"Glu-";

			else if ((codon.charAt(0)=='G' && codon.charAt(1)=='G'))
				protein=protein+"Gly-";

			else if (codon.equals("GAA")||codon.equals("GAG"))
				protein=protein+"Glu-";

			else throw new IllegalArgumentException("One or more bases are incorrect");
		}
		System.out.println(protein.substring(0, protein.length()-1));
	}

	
	/**
	 * @return the index of the start codon
	 */
	public static int start(String str){
		int i=0;
		while (i<str.length()-2){
			if(str.charAt(i)=='A' && str.charAt(i+1)=='U' && str.charAt(i+2)=='G')
				return i;
			else
				i++;
		}
		throw new IllegalArgumentException("\n START codon (AUG) is missing");
	}

	/**
	 * @return the index of the end codon
	 */
	public static int end(String str){ 
		int i=start(str);
		while (i<str.length()-2){
			if((str.charAt(i)=='U' && str.charAt(i+1)=='A' && str.charAt(i+2)=='A')||
					(str.charAt(i)=='U' && str.charAt(i+1)=='A' && str.charAt(i+2)=='G')||
					(str.charAt(i)=='U' && str.charAt(i+1)=='G' && str.charAt(i+2)=='A'))
				return i;
			else
				i=i+3;
		}
		throw new IllegalArgumentException("\n END codon (UAA/UAG/UGA) is missing");
	}
}
