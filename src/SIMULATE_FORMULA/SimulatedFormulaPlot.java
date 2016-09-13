package SIMULATE_FORMULA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;

import org.openscience.cdk.formula.IsotopeContainer;

import MISC.ToolBox;

public class SimulatedFormulaPlot {

	public static void execute(String[] args) {
		
		try {
			String filePath = args[0];
			File folder = new File(filePath);
			
			FileWriter fwriter = new FileWriter(args[1]); 
			BufferedWriter out = new BufferedWriter(fwriter);
			out.write("Mass\tCountC\tCountN\tCountCN\tCountS\tCountP\tCountPassHydrogenRule\tCountFailedHydrogenRule\n");
			int count = 0;
			File[] listOfFiles = folder.listFiles();
			for (int i = 0; i < listOfFiles.length; i++) {
				if (listOfFiles[i].isFile()) {
					int countC = 0;
					int countN = 0;
					int countS = 0;
					int countP = 0;
					int countCN = 0;	
					int countPassHydrogenRule = 0;
					int countFailedHydrogenRule = 0;
					
					FileInputStream fstream = new FileInputStream(listOfFiles[i].getPath());
					DataInputStream din = new DataInputStream(fstream);
					BufferedReader in = new BufferedReader(new InputStreamReader(din));
					while (in.ready()) {
						String str = in.readLine();
						String[] split = str.split(":");
						String formula = split[0];
						String mass = split[1];
						if (!(formula.contains("Br") || formula.contains("F") || formula.contains("Cl"))) {
							if (formula.contains("C")) {
								countC++;
							}
							if (formula.contains("N")) {
								countN++;
							}
							if (formula.contains("S")) {
								countS++;
							}
							if (formula.contains("P")) {
								countP++;
							}
							if (formula.contains("C") && formula.contains("N")) {
								countCN++;
							}
							if (ToolBox.check_hydrogen_rule(formula)) {
								countPassHydrogenRule++;
							} else {
								countFailedHydrogenRule++;
							}
						}
					}
					in.close();
					out.write(listOfFiles[i].getName().replaceAll("Formulas", "") + "\t" + countC + "\t" + countN + "\t" + countCN + "\t" + countP + "\t" + countS + "\t" + countPassHydrogenRule + "\t" + countFailedHydrogenRule + "\n");
				}
			}
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
