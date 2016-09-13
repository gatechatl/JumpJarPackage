package MassDatabaseGeneration;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import MISC.ToolBox;

/**
 * After the mass query database is generated
 * To assist the querying, we also need to add other database information.
 * @author tshaw
 *
 */
public class AppendDatabaseInfo {
 
	
	public static void execute(String[] args) {
		
		try {
			HashMap database_formula = new HashMap();
			
			String referenceFile = args[0];
			int massColIndex = new Integer(args[1]);
			String databasePath = args[2];
			String databaseName = args[3];
			String test_flag_str = args[4];
			boolean test_flag = false;
			if (test_flag_str.toUpperCase().equals("TRUE")) {			
				test_flag = true;			
			}
			int id = 0;
			FileInputStream fstream = new FileInputStream(referenceFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String formula = in.readLine().trim();
				String[] split = formula.split("\t");
				if (split.length > massColIndex) {
					if (ToolBox.check_formula_valid_element(split[massColIndex].split(":")[0])) {
						
						double mass = ToolBox.getMonoisotopicMass(split[massColIndex].split(":")[0]);
						database_formula.put(ToolBox.HillSystemOrder(split[massColIndex].split(":")[0]), mass);
					}
				}
			}
			in.close();
			
			Iterator itr = database_formula.keySet().iterator();
			while (itr.hasNext()) {
				String formula = (String)itr.next();
				double mass = (Double)database_formula.get(formula);
				int index = new Double(mass * 100).intValue() + 1;
				
				String formula_file = databasePath + "/Formula" + index + ".txt";
				File file = new File(formula_file);
				if (file.exists()) {
					LinkedList list = new LinkedList();
					id++;
					FileInputStream formula_fstream = new FileInputStream(formula_file);
					DataInputStream formula_din = new DataInputStream(formula_fstream);
					BufferedReader formula_in = new BufferedReader(new InputStreamReader(formula_din));
					while (formula_in.ready()) {
						String original_formula_mass = formula_in.readLine().trim();
						String[] split = original_formula_mass.split("\t")[0].split(":");
						String original_formula = ToolBox.HillSystemOrder(split[0]);
						String original_mass = split[1];
						if (original_formula.equals(formula)) {
							String database_info = "";
							String[] split2 = original_formula_mass.split("\t");
							if (split2.length == 2) {
								database_info = split2[1];
								String[] separate_existing_db = database_info.split(",");
								boolean existing = false;
								for (String db: separate_existing_db) {
									if (db.equals(databaseName)) {
										existing = true;
									}
								}
								if (!existing) {
									database_info += databaseName + ",";
								}
							} else {
								database_info += databaseName + ",";
							}
							list.add(split2[0] + "\t" + database_info);
						} else {
							list.add(original_formula_mass);
						}
						if (test_flag) {
							while (formula_in.ready()) {
								formula_in.readLine();
							}
						}
					}
					formula_in.close();
					
					if (!test_flag) {
						String outputFile = formula_file;
						FileWriter fwriter = new FileWriter(outputFile);
						BufferedWriter out = new BufferedWriter(fwriter);					
						Iterator itr2 = list.iterator();
						while (itr2.hasNext()) {
							String str = (String)itr2.next();
							out.write(str + "\n");
						}
						out.close();
					} else {
						Iterator itr2 = list.iterator();
						while (itr2.hasNext()) {
							String str = (String)itr2.next();
							System.out.println(str);
						}
					}
				} else {
					System.out.println("Missing: " + formula_file + "\t" + id);
				}
			
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
