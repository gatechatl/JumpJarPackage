package PUBCHEM_STRUCTUREDATABASE;

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

import CDKFunction.CalculateMonoisotope;
import METABOLOMICS.AIM.QueryStructureDatabase;
import MISC.ToolBox;

/**
 * Go through each formula and add annotation to each structure
 * @author tshaw
 *
 */
public class PUBCHEMStructureAddAnnotation {

	public static String parameter_info() {
		return "[hmdb_file] [kegg_file] [ymdb_file] [database]";
	}
	
	public static void main(String[] args) {
		
		String str = "InChI=1S/C17H16O5/c1-10(2)6-8-21-17-15-12(7-9-20-15)14(19-3)11-4-5-13(18)22-16(11)17/h4-7,9H,8H2,1-3H3";
		if (str.contains("InChI=")) {
			System.out.println("FOund");
			String elementOnly = QueryStructureDatabase.elementOnly("C21H34O5");
			boolean foo = ToolBox.check_formula_valid_element("C21H34O5");
			double monoisotopicMass = CalculateMonoisotope.getMonoisotopicMassCDK("C21H34O5");
			System.out.println(elementOnly);
			System.out.println(monoisotopicMass);
			System.out.println(foo);
		}
	}
	public static void execute(String[] args) {
		
		try {
			
			String formula_file = args[0];
			String kegg_file = args[1];
			String ymdb_file = args[2];
			String database = args[3];
			
			HashMap map = readHMDB(formula_file);
			HashMap kegg_map = readKEGG(kegg_file);
			HashMap ymdb_map = readYMDB(ymdb_file);
			
			//String datatype = args[2];
			
			System.out.println("YMDB Map: " + ymdb_map.size());
			Iterator itr = ymdb_map.keySet().iterator();
			while (itr.hasNext()) {
				String ymdb_inchikey = (String)itr.next();
				YMDB_DATA data = (YMDB_DATA)ymdb_map.get(ymdb_inchikey);
				String elementOnly = QueryStructureDatabase.elementOnly(data.FORMULA);
				//System.out.println(data.FORMULA);
				
				double monoisotopicMass = CalculateMonoisotope.getMonoisotopicMassCDK(data.FORMULA);
				int mass = new Double(monoisotopicMass).intValue();
				//System.out.println(mass + "\t" + elementOnly);
				
				HashMap list = new HashMap();
				
				String fileName = database + "/" + mass + "/" + elementOnly + ".txt";
				File f = new File(database + "/" + mass + "/" + elementOnly + ".txt");
								
				System.out.println(fileName);
				//System.out.println("Mass is this: " + mass);
				if (f.isFile()) {
					System.out.println("Search: " + ymdb_inchikey + "\t" + data.FORMULA + "\t" + mass);
					LinkedList original = new LinkedList();
					FileInputStream fstream = new FileInputStream(fileName);
					DataInputStream din = new DataInputStream(fstream);
					BufferedReader in = new BufferedReader(new InputStreamReader(din));
					while (in.ready()) {
						String str = in.readLine();
						original.add(str);
					}
					in.close();
				
					FileWriter fwriter = new FileWriter(database + "/" + mass + "/" + elementOnly + ".txt");
					BufferedWriter out = new BufferedWriter(fwriter);
					
					boolean hit = false;
					Iterator itr2 = original.iterator();
					while (itr2.hasNext()) {
						String str = (String)itr2.next();
						String[] split = str.split("\t");
						if (split.length > 4) {
							String database_formula = split[1];
							String inchikey = split[2];						
							String inchi = split[3];	
							String smile = split[4];
							String flag = "";
							String name = "";
							String clean_hmdbinchikey = ymdb_inchikey.replaceAll("InChIKey=", "");
							
							if (inchikey.equals(clean_hmdbinchikey)) {
								//System.out.println("found");
								if (split.length >= 9) {
									name = split[7];
									flag = split[8];
									if (name.contains("InChI=")) {
										name = "\"" + data.NAME + "\"";
									}
									
									if (!flag.contains("YMDB")) {
										flag += ",YMDB";
									}
									if (!name.equals("")) {
										if (!name.contains(data.NAME)){ 
											name += ",\"" + data.NAME + "\"";
										}
									}
								} else {
									flag = "YMDB";
									name = "\"" + data.NAME + "\"";
								}
								String found = split[0] + "\t" + split[1] + "\t" + split[2] + "\t" + split[3] + "\t" + split[4] + "\t" + split[5] + "\t" + split[6] + "\t" + name + "\t" + flag;
								hit = true;
								System.out.println("YMDB_Found:");
								System.out.println(found);
								out.write(found + "\n");
							} else {
								//System.out.println("YMDB_NoFound");
								out.write(str + "\n");
							}
							
							
							
						}
					} // end iterator

					
					if (!hit) {
							System.out.println("Not Found: ");
							String notfound_str = data.ID + "\t" + data.FORMULA + "\t" + data.INCHIKEY.replaceAll("InChIKey=", "") + "\t" + data.INCHI.replaceAll("InChI=", "") + "\t" + data.SMILE + "\t" + data.IUPAC + "\t" + data.MASS + "\t" + data.NAME + "\tYMDB"; 
							System.out.println(notfound_str);
							out.write(notfound_str + "\n");
					}
					
					out.close();
				} //
			} // end ymdb
			
			System.out.println("Kegg Map: " + kegg_map.size());
			itr = kegg_map.keySet().iterator();
			while (itr.hasNext()) {
				String kegg_inchikey = (String)itr.next();
				KEGG_DATA data = (KEGG_DATA)kegg_map.get(kegg_inchikey);
				String elementOnly = QueryStructureDatabase.elementOnly(data.FORMULA);
				//System.out.println(data.FORMULA);
				double monoisotopicMass = CalculateMonoisotope.getMonoisotopicMassCDK(data.FORMULA);
				int mass = new Double(monoisotopicMass).intValue();
				//System.out.println(mass + "\t" + elementOnly);
				
				HashMap list = new HashMap();
				
				String fileName = database + "/" + mass + "/" + elementOnly + ".txt";
				File f = new File(database + "/" + mass + "/" + elementOnly + ".txt");
				
				
				System.out.println("Kegg\t" + kegg_inchikey + "\t" + fileName);
				//System.out.println("Mass is this: " + mass);
				if (f.isFile()) {
					System.out.println("Search: " + kegg_inchikey + "\t" + data.FORMULA + "\t" + mass);
					LinkedList original = new LinkedList();
					FileInputStream fstream = new FileInputStream(fileName);
					DataInputStream din = new DataInputStream(fstream);
					BufferedReader in = new BufferedReader(new InputStreamReader(din));
					while (in.ready()) {
						String str = in.readLine();
						original.add(str);
					}
					in.close();
				
					FileWriter fwriter = new FileWriter(database + "/" + mass + "/" + elementOnly + ".txt");
					BufferedWriter out = new BufferedWriter(fwriter);
					
					boolean hit = false;
					Iterator itr2 = original.iterator();
					while (itr2.hasNext()) {
						String str = (String)itr2.next();
						String[] split = str.split("\t");
						if (split.length >= 7) {
							String database_formula = split[1];
							String inchikey = split[2];						
							String inchi = split[3];	
							String smile = split[4];
							String flag = "";
							String name = "";
							String clean_hmdbinchikey = kegg_inchikey.replaceAll("InChIKey=", "");
							
							if (inchikey.equals(clean_hmdbinchikey)) {
								//System.out.println("found");
								if (split.length >= 9) {
									name = split[7];
									if (name.contains("InChI=")) {
										name = "";
										name = "\"" + data.NAME + "\"";
									}
									//name = cleanNAME(name);
									flag = split[8];
									
									if (!flag.contains("KEGG")) {
										flag += ",KEGG";
									}
									if (!name.equals("")) {
										if (!name.contains(data.NAME)){ 
											name += ",\"" + data.NAME + "\"";
										}
									}
								} else {
									flag = "KEGG";
									name = "\"" + data.NAME + "\"";
								}
								String found = split[0] + "\t" + split[1] + "\t" + split[2] + "\t" + split[3] + "\t" + split[4] + "\t" + split[5] + "\t" + split[6] + "\t" + name + "\t" + flag;
								hit = true;
								System.out.println("Found:");
								System.out.println(found);
								out.write(found + "\n");
							} else {
								out.write(str + "\n");
							}
							
							
							
						}
					} // end iterator

					
					if (!hit) {
							System.out.println("Not Found: ");
							String notfound_str = data.ID + "\t" + data.FORMULA + "\t" + data.INCHIKEY.replaceAll("InChIKey=", "") + "\t" + data.INCHI.replaceAll("InChI=", "") + "\t" + data.SMILE + "\t" + data.IUPAC + "\t" + data.MASS + "\t" + data.NAME + "\tKEGG"; 
							System.out.println(notfound_str);
							out.write(notfound_str + "\n");
					}
					
					out.close();
				} //
			} // end kegg
			
			
			
			itr = map.keySet().iterator();
			while (itr.hasNext()) {
				String hmdb_inchikey = (String)itr.next();
				HMDB_DATA data = (HMDB_DATA)map.get(hmdb_inchikey);
				String elementOnly = QueryStructureDatabase.elementOnly(data.FORMULA);
				//System.out.println(data.FORMULA);
				double monoisotopicMass = CalculateMonoisotope.getMonoisotopicMassCDK(data.FORMULA);
				int mass = new Double(monoisotopicMass).intValue();
				//System.out.println(mass + "\t" + elementOnly);
				
				HashMap list = new HashMap();
				
				String fileName = database + "/" + mass + "/" + elementOnly + ".txt";
				File f = new File(database + "/" + mass + "/" + elementOnly + ".txt");
				
				
				//System.out.println(fileName);
				//System.out.println("Mass is this: " + mass);
				if (f.isFile()) {
					System.out.println("Search: " + hmdb_inchikey + "\t" + data.FORMULA + "\t" + mass);
					LinkedList original = new LinkedList();
					FileInputStream fstream = new FileInputStream(fileName);
					DataInputStream din = new DataInputStream(fstream);
					BufferedReader in = new BufferedReader(new InputStreamReader(din));
					while (in.ready()) {
						String str = in.readLine();
						original.add(str);
					}
					in.close();
				
					FileWriter fwriter = new FileWriter(database + "/" + mass + "/" + elementOnly + ".txt");
					BufferedWriter out = new BufferedWriter(fwriter);
					
					boolean hit = false;
					Iterator itr2 = original.iterator();
					while (itr2.hasNext()) {
						String str = (String)itr2.next();
						String[] split = str.split("\t");
						if (split.length > 4) {
							String database_formula = split[1];
							String inchikey = split[2];						
							String inchi = split[3];	
							String smile = split[4];
							String flag = "";
							String name = "";
							String clean_hmdbinchikey = hmdb_inchikey.replaceAll("InChIKey=", "");
							
							if (inchikey.equals(clean_hmdbinchikey)) {
								//System.out.println("found");
								if (split.length >= 9) {
									name = split[7];
									//name = cleanNAME(name);
									flag = split[8];
									
									if (!flag.contains("HMDB")) {
										flag += ",HMDB";
									}
									if (!name.equals("")) {
										if (!name.contains(data.NAME)){ 
											name += ",\"" + data.NAME + "\"";
										}
									}
								} else {
									flag = "HMDB";
									name = "\"" + data.NAME + "\"";
								}
								String found = split[0] + "\t" + split[1] + "\t" + split[2] + "\t" + split[3] + "\t" + split[4] + "\t" + split[5] + "\t" + split[6] + "\t" + name + "\t" + flag;
								hit = true;
								System.out.println("Found:");
								System.out.println(found);
								out.write(found + "\n");
							} else {
								out.write(str + "\n");
							}
						}
					} // end iterator
					
					if (!hit) {
							System.out.println("Not Found: ");
							String notfound_str = data.ID + "\t" + data.FORMULA + "\t" + data.INCHIKEY.replaceAll("InChIKey=", "") + "\t" + data.INCHI.replaceAll("InChI=", "") + "\t" + data.SMILE + "\t" + data.IUPAC + "\t" + data.MASS + "\t" + data.NAME + "\tHMDB"; 
							System.out.println(notfound_str);
							out.write(notfound_str + "\n");
					}
					out.close();
				} //
			} // end hmdb
									
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static String cleanNAME(String name) {
		String result = "";
		String[] split = name.split(",");
		for (String str: split) {
			if (!str.contains("InChI=")) {
				if (result.equals("")) {
					result = str;
				} else {
					result += "," + str;
				}
			}
		}
		return result;
	}
	public static String cleanSMILE(String hmdb_smile) {
		return hmdb_smile.replaceAll("\\\\", "").replaceAll("\\/", "");
	}
	public static HashMap readYMDB(String inputFile) {
		HashMap map = new HashMap();
		try {
			FileInputStream fstream2 = new FileInputStream(inputFile);
			DataInputStream din2 = new DataInputStream(fstream2);
			BufferedReader in2 = new BufferedReader(new InputStreamReader(din2));
			while (in2.ready()) {
				String str2 = in2.readLine();
				String[] split = str2.split("\t");
				YMDB_DATA data = new YMDB_DATA();
				System.out.println(split[0] + "!\t" + str2);
				if (split.length > 6) {
					data.ID = split[0];
					data.FORMULA = split[1];
					data.INCHIKEY = split[2];
					data.INCHI = split[3];
					data.SMILE = split[4];
					data.IUPAC = split[5];
					data.MASS = split[6];
					data.NAME = split[5];
					data.ROW = str2;
					if (split.length >= 7) {
						
						if (ToolBox.check_formula_valid_element(data.FORMULA)) {
							map.put(data.INCHIKEY, data);
							//System.out.println(data.ID + "\t" + data.FORMULA);
						}
					}
				}
			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return map;
	}
	public static HashMap readKEGG(String inputFile) {
		HashMap map = new HashMap();
		try {
			FileInputStream fstream2 = new FileInputStream(inputFile);
			DataInputStream din2 = new DataInputStream(fstream2);
			BufferedReader in2 = new BufferedReader(new InputStreamReader(din2));
			while (in2.ready()) {
				String str2 = in2.readLine();
				String[] split = str2.split("\t");
				KEGG_DATA data = new KEGG_DATA();
				
				if (split.length > 6) {
					data.ID = split[0];
					data.FORMULA = split[1];
					data.INCHIKEY = split[2];
					data.INCHI = split[3];
					data.SMILE = split[4];
					data.IUPAC = split[5];
					data.MASS = split[6];
					data.NAME = split[5];
					data.ROW = str2;
					if (split.length >= 7) {
						
						if (ToolBox.check_formula_valid_element(data.FORMULA)) {
							map.put(data.INCHIKEY, data);
							System.out.println(split[0] + "!\t" + str2);
							//System.out.println(data.ID + "\t" + data.FORMULA);
						}
					}
				}
			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return map;
	}
	public static HashMap readHMDB(String inputFile) {
		HashMap map = new HashMap();
		try {
			FileInputStream fstream2 = new FileInputStream(inputFile);
			DataInputStream din2 = new DataInputStream(fstream2);
			BufferedReader in2 = new BufferedReader(new InputStreamReader(din2));
			while (in2.ready()) {
				String str2 = in2.readLine();
				String[] split = str2.split("\t");
				HMDB_DATA data = new HMDB_DATA();
				data.ID = split[0];
				data.FORMULA = split[1];
				data.INCHIKEY = split[2];
				data.INCHI = split[3];
				data.SMILE = split[4];
				data.IUPAC = split[5];
				data.MASS = split[6];
				data.NAME = split[split.length - 1];
				data.ROW = str2;
				if (split.length >= 8) {
					
					if (ToolBox.check_formula_valid_element(data.FORMULA)) {
						map.put(data.INCHIKEY, data);
						//System.out.println(data.ID + "\t" + data.FORMULA);
					}
				}
				
			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return map;
	}
	public static class YMDB_DATA {
		public String ID = "";
		public String FORMULA = "";
		public String INCHIKEY = "";
		public String INCHI = "";
		public String SMILE = "";
		public String IUPAC = "";
		public String MASS = "";
		public String NAME = "";
		public String ROW = "";
	}
	public static class KEGG_DATA {
		public String ID = "";
		public String FORMULA = "";
		public String INCHIKEY = "";
		public String INCHI = "";
		public String SMILE = "";
		public String IUPAC = "";
		public String MASS = "";
		public String NAME = "";
		public String ROW = "";
	}
	public static class HMDB_DATA {
		public String ID = "";
		public String FORMULA = "";
		public String INCHIKEY = "";
		public String INCHI = "";
		public String SMILE = "";
		public String IUPAC = "";
		public String MASS = "";
		public String NAME = "";
		public String ROW = "";
	}
}
