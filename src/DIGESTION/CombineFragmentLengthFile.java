package DIGESTION;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class CombineFragmentLengthFile {

	public static String parameter_info() {
		return "[inputFile1] [inputFile2] [Alias1] [Alias2]";
	}
	public static void execute(String[] args) {
		
		try {
			
			String inputFile1 = args[0];
			String inputFile2 = args[1];
			String alias1 = args[2];
			String alias2 = args[3];
			
			System.out.println("Type\tLength");
			FileInputStream fstream = new FileInputStream(inputFile1);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				System.out.println(alias1 + "\t" + str);
			}
			in.close();
			
			fstream = new FileInputStream(inputFile2);
			din = new DataInputStream(fstream);
			in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				System.out.println(alias2 + "\t" + str);
			}
			in.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
