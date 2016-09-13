package MassDatabaseGeneration;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.LinkedList;

/**
 * Go through the database and filter out formulas
 * @author tshaw
 *
 */
public class FilterFClBr {

	public static String parameter_info() {
		return "[OldDatabasePath] [NewDatabasePath]";
	}
	public static void execute(String[] args) {
		try {
			String databasePath = args[0];
			String newDatabasePath = args[1];
			
			File folder = new File(newDatabasePath);
			if (!folder.exists()) {
				if (!folder.isDirectory()) {
					folder.mkdir();
				}
			}
			
			for (int i = 0; i < 150100; i++) {
				String index = i + "";
				
				
				String formula_file = databasePath + "/Formula" + index + ".txt";
				File file = new File(formula_file);
				if (file.exists()) {
					String new_formula_file = newDatabasePath + "/Formula" + index + ".txt";
					FileWriter fwriter = new FileWriter(new_formula_file); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\All_Formulas.txt");
					BufferedWriter out = new BufferedWriter(fwriter);
					
					FileInputStream formula_fstream = new FileInputStream(formula_file);
					DataInputStream formula_din = new DataInputStream(formula_fstream);
					BufferedReader formula_in = new BufferedReader(new InputStreamReader(formula_din));
					while (formula_in.ready()) {
						String str = formula_in.readLine();
						String[] split = str.split("\t");
						if (!(split[0].contains("Cl") || split[0].contains("Br") || split[0].contains("F"))) {
							out.write(str + "\n");
						}						
					}
					formula_in.close();
					out.close();
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}	
	}
	
	
}
