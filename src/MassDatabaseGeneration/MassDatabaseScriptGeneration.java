package MassDatabaseGeneration;

import java.io.BufferedWriter;
import java.io.FileWriter;

public class MassDatabaseScriptGeneration {

	public static void execute(String[] args) {
		try {
			
			String output_mass_prefix = args[0];
			String output_shell_prefix = args[1];
			//String output_final_shell = args[2];
			int num = new Integer(args[2]);
			String formulaType = args[3];
			//
			int n = num; 
			
			int count = 0;
			
			
			FileWriter fwriter2 = new FileWriter(output_shell_prefix); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\All_Formulas.txt");
			BufferedWriter out2 = new BufferedWriter(fwriter2);
			//FileWriter fwriter3 = new FileWriter(output_final_shell); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\All_Formulas.txt");
			//BufferedWriter out3 = new BufferedWriter(fwriter3);
			
			for (int j = 1; j <= n; j++) {
				FileWriter fwriter = new FileWriter(output_mass_prefix + "_" + j); 
				BufferedWriter out = new BufferedWriter(fwriter);

				//out2.write("java -cp .:lib/cdk-1.4.19.jar MISSILE.GetAllMetaboliteFormulas Mass_" + j + " 0.1 FormulaNumMass_" + j + ".txt Formulas_" + j + ".txt false\n");
				out2.write("metdbgen Mass_" + j + " 0.01 FormulaNumMass_" + j + ".txt Formulas_" + j + ".txt " + formulaType + "\n");
				
				//out3.write("bsub -P " + j + " -eo error/simple" + j + ".err.txt -oo error/simple" + j + ".out.txt \"sh " + args[1] + "_" + j + "\"\n");
				for (double i = 0.01 * j; i <= 1500.0001; i = i + 0.01 * n) {
					//System.out.println(i);
					
					out.write(i + "\n");
					count++;
				}
				out.close();
			}
			out2.close();
			//out3.close();
			System.out.println(count);
			//
		} catch (Exception e) {
			e.printStackTrace();
			
		}
	}
}


