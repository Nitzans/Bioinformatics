import java.util.Scanner;

public class RNAtranscription {
	
	/**
	 * @Description:
	 * Received from to user DNA strand and translate it to suitable RNA strand which will be created from
	 */
	public static void main(String[] args){
		String DNA;
		String RNA = "";
		System.out.println("Enter DNA strand: ");
		Scanner sc = new Scanner(System.in);
		DNA = sc.next();
		
		for(int i=0;i<DNA.length();i++){
			if (DNA.charAt(i)=='A')
				RNA= RNA+"U";
			else if (DNA.charAt(i)=='T')
				RNA= RNA+"A";
			else if (DNA.charAt(i)=='C')
				RNA= RNA+"G";
			else if (DNA.charAt(i)=='G')
				RNA= RNA+"C";
			else throw new IllegalStateException("\n Please enter valid base A/T/C/G and try again");
		}
		System.out.println(RNA);
	}
}
