package MassDatabaseGeneration;

import java.io.BufferedWriter;
import java.io.FileWriter;

public class CreateIndexFormulaScript {

	public static void execute(String[] args) {
		
		try {
			
			String inputFile = args[0];
			String outputFolder = args[1];
			int min = new Integer(args[2]);
			int max = new Integer(args[3]);
			int freq = new Integer(args[4]);
			String outputFile = args[5];
			FileWriter fwriter = new FileWriter(outputFile);
			BufferedWriter out = new BufferedWriter(fwriter);
			
			int diff = max - min;
			for (int i = 0; i < freq; i++) {
				int low = min + (diff / freq) * i;
				int high = min + (diff / freq) * (i + 1);
				
				System.out.println("jumpjar -CreateIndexFormula " + inputFile + " " + outputFolder + " " + low + " " + high);
				out.write("jumpjar -CreateIndexFormula " + inputFile + " " + outputFolder + " " + low + " " + high + "\n");
			}
			out.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
