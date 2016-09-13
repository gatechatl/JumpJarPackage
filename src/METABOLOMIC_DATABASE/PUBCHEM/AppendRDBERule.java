package METABOLOMIC_DATABASE.PUBCHEM;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;

import METABOLOMIC_DATABASE.PUBCHEM.DownloadPubchemXML.PUBCHEM;
import MISC.RDBERule;
import MISC.ToolBox;

public class AppendRDBERule {

	
	public static String parameter_info() {
		return "[xmlFile] [outputFile]";
	}
	public static void execute(String[] args) {
		
		try {
			
			//String inputFile = args[0];
			String inputPath = args[0];
			String outputPath = args[1];
			
			File f = new File(inputPath);
			for (File file: f.listFiles()) {
				
				String fileName = file.getPath();
				System.out.println("Reading file: " + fileName);
				if (!fileName.contains(".out")) {
					String outputFile = outputPath + "/" + file.getName();
					File f2 = new File(outputFile);
					if (f2.exists()) {
						System.out.println("File already exist");
						System.exit(0);
					}
					FileWriter fwriter = new FileWriter(outputFile);
		            BufferedWriter out = new BufferedWriter(fwriter);
					
					FileInputStream fstream = new FileInputStream(fileName);
					DataInputStream din = new DataInputStream(fstream);
					BufferedReader in = new BufferedReader(new InputStreamReader(din));
					while (in.ready()) {
						String str = in.readLine();
						String[] split = str.split("\t");
						String formula = split[1];
						boolean rdbe = false;
						formula = formula.replaceAll("\\-", "");
						formula = formula.replaceAll("\\+", "");
						//System.out.print(formula);
						formula = ToolBox.standardize_name(formula);
						//System.out.println("\t" + formula);
						boolean check_formula_valid = ToolBox.check_formula_valid_element_corrected(formula);
						boolean check_hydrogen_rule = ToolBox.check_hydrogen_rule(formula);
						RDBERule rule = new RDBERule();

						//System.out.println(formula);
						if (rule.validateRBDE(ToolBox.standardize_name(formula))) {
							rdbe = true;
						}
						out.write(split[0] + "\t" + split[1] + "\t" + split[2] + "\t" + split[3] + "\t" + split[4] + "\t" + split[5] + "\t" + split[6] + "\t" + split[7] + "\t" + split[8] + "\t" + split[9] + "\t" + check_formula_valid + "\t" + check_hydrogen_rule + "\t" + rdbe + "\n");
						//out.write(str + "\t" + rbde + "\n");
					}
					in.close();		           		            
					out.close();
					
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
