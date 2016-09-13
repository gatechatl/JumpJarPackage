package MassDatabaseGeneration;

import java.io.BufferedWriter;
import java.io.FileWriter;

public class CreateSimulatedRangeOld {

	public static void main(String[] args) {
		try {
			
			String output_prefix = args[0];
			
			//
			int n = 60;
			
			int count = 0;
			
			FileWriter fwriter3 = new FileWriter(args[2]); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\All_Formulas.txt");
			BufferedWriter out3 = new BufferedWriter(fwriter3);
			
			for (int j = 1; j <= n; j++) {
				FileWriter fwriter = new FileWriter(output_prefix + "_" + j); 
				BufferedWriter out = new BufferedWriter(fwriter);
				
				FileWriter fwriter2 = new FileWriter(args[1] + "_" + j); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\All_Formulas.txt");
				BufferedWriter out2 = new BufferedWriter(fwriter2);
				out2.write("java -cp .:lib/cdk-1.4.19.jar MISSILE.GetAllMetaboliteFormulas Mass_" + j + " 0.1 FormulaNumMass_" + j + ".txt Formulas_" + j + ".txt false\n");
				out2.close();
				out3.write("bsub -P " + j + " -eo error/simple" + j + ".err.txt -oo error/simple" + j + ".out.txt \"sh " + args[1] + "_" + j + "\"\n");
				for (double i = 0.01 * j; i <= 1500.0001; i = i + 0.01 * n) {
					//System.out.println(i);
					
					out.write(i + "\n");
					count++;
				}
				out.close();
			}
			out3.close();
			System.out.println(count);
			//
		} catch (Exception e) {
			e.printStackTrace();
			
		}
	}
}

