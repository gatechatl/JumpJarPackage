package PUBCHEM_STRUCTUREDATABASE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;

import MISC.ToolBox;

/**
 * 
 * @author tshaw
 *
 */
public class AppendPubChemTable {

	public static String description() {
		return "append hmdb ymdb to pubchem";
	}
	public static String parameter_info() {
		return "[hmdb_file] [ymdb_file] [pubchem_file] [outputFile]";
	}
	public static void execute(String[] args) {
		
		try {
			
			HashMap hmdb_map = new HashMap();
			HashMap hmdb_hit_map = new HashMap();
			HashMap ymdb_map = new HashMap();
			HashMap ymdb_hit_map = new HashMap();
			String hmdb_fileName = args[0];
			String ymdb_fileName = args[1];
			String pubchem_fileName = args[2];
			String outputFile = args[3];
			String problemFile = args[4];
			
        	FileWriter fwriter = new FileWriter(outputFile);
            BufferedWriter out = new BufferedWriter(fwriter);
            
            FileWriter fwriter2 = new FileWriter(problemFile);
            BufferedWriter out2 = new BufferedWriter(fwriter2);
            
			FileInputStream fstream = new FileInputStream(hmdb_fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (split.length > 5) {
					String id = split[0];
					String formula = split[1];
					String inchikey = split[2].replaceAll("InChIKey=", "");
					String inchi = split[3].replaceAll("InChI=", "");
					String smile = split[4];
					String iupac = split[5];
					double mass = ToolBox.getMonoisotopicMass(formula);
					String regularName = split[7];
					hmdb_map.put(inchikey, id + "\t" + formula + "\t" + inchikey + "\t" + inchi + "\t" + smile + "\t" + iupac + "\t" + mass + "\t" + regularName);
				}
			}
			in.close();
			
			fstream = new FileInputStream(ymdb_fileName);
			din = new DataInputStream(fstream);
			in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				String id = split[0];
				if (split.length > 5) {
					String formula = split[1];
					String inchikey = split[2].replaceAll("InChIKey=", "");
					String inchi = split[3].replaceAll("InChI=", "");
					String smile = split[4];
					String iupac = split[5];
					double mass = ToolBox.getMonoisotopicMass(formula);
					String regularName = "NA";
				
					ymdb_map.put(inchikey, id + "\t" + formula + "\t" + inchikey + "\t" + inchi + "\t" + smile + "\t" + iupac + "\t" + mass + "\t" + regularName);
				}
			}
			in.close();
			
			double count = 0;
			fstream = new FileInputStream(pubchem_fileName);
			din = new DataInputStream(fstream);
			in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();				
				String[] split = str.split("\t");
				if (split.length > 5) {
					String pubchem_id = split[0];
					String formula = split[1];
					String inchikey = split[2].replaceAll("InChIKey=", "");
					String inchi = split[3].replaceAll("InChI=", "");
					String smile = split[4];
					String iupac = split[5];
					double mass = ToolBox.getMonoisotopicMass(formula);
					String type = "PUBCHEM";
					String general_name = "NA";
					if (hmdb_map.containsKey(inchikey)) {
						type += ",HMDB";
						String line = (String)hmdb_map.get(inchikey);
						hmdb_hit_map.put(inchikey, line);
						String[] split2 = line.split("\t");
						general_name = split2[7];
					}
					if (ymdb_map.containsKey(inchikey)) {
						String line = (String)ymdb_map.get(inchikey);
						ymdb_hit_map.put(inchikey, line);
						type += ",YMDB";
						/*if (type.equals("PUBCHEM")) {
							type += ",YMDB";
						}*/
					}
					out.write(pubchem_id + "\t" + formula + "\t" + inchikey + "\t" + inchi + "\t" + smile + "\t" + iupac + "\t" + mass + "\t" + general_name + "\t" + type + "\n");
					count++;
					if (count % 100000 == 0) {
						out.flush();
					}
				} else {
					out2.write(str + "\n");
				}
			}
			in.close();
			
			
			Iterator itr = hmdb_map.keySet().iterator();
			while (itr.hasNext()) {
				String inchikey = (String)itr.next();
				String general_name = "NA";
				String type = "";
				if (!hmdb_hit_map.containsKey(inchikey)) {
					String line = (String)hmdb_map.get(inchikey);
					String[] split = line.split("\t");
					general_name = split[7];
					type = "HMDB";
					if (ymdb_map.containsKey(inchikey)) {
						type += "YMDB";
					} 
					String hmdb_id = split[0];
					String formula = split[1];					
					String inchi = split[3].replaceAll("InChI=", "");
					String smile = split[4];
					String iupac = split[5];
					double mass = ToolBox.getMonoisotopicMass(formula);
					out.write(hmdb_id + "\t" + formula + "\t" + inchikey + "\t" + inchi + "\t" + smile + "\t" + iupac + "\t" + mass + "\t" + general_name + "\t" + type + "\n");
				} 																
			}
			
			
			itr = ymdb_map.keySet().iterator();
			while (itr.hasNext()) {
				String inchikey = (String)itr.next();
				String general_name = "NA";
				String type = "";
				if (!hmdb_map.containsKey(inchikey) && !ymdb_hit_map.containsKey(inchikey)) {
					String line = (String)ymdb_map.get(inchikey);
					String[] split = line.split("\t");
					general_name = split[7];
					type = "YMDB";
					
					String ymdb_id = split[0];
					String formula = split[1];					
					String inchi = split[3].replaceAll("InChI=", "");
					String smile = split[4];
					String iupac = split[5];
					double mass = ToolBox.getMonoisotopicMass(formula);
					out.write(ymdb_id + "\t" + formula + "\t" + inchikey + "\t" + inchi + "\t" + smile + "\t" + iupac + "\t" + mass + "\t" + general_name + "\t" + type + "\n");
				} 																
			}
			
			out.close();
			out2.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
