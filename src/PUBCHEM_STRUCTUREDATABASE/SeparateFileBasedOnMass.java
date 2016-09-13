package PUBCHEM_STRUCTUREDATABASE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import MISC.ToolBox;

public class SeparateFileBasedOnMass {

	public static String parameter_info() {
		return "[inputFile] [size_mass_range] [outputFile] [outputScript] [jumpjarPath] [outputFolder]";
	}
	public static void execute(String[] args) {
		
		try {
			
			String inputFile = args[0];			
			int num = new Integer(args[1]);
			String outputFile = args[2];
			String outputScript = args[3];
			String jumpjarPath = args[4];
			String outputFolder = args[5];
			HashMap map = new HashMap();
			LinkedList list = new LinkedList();
			int idx = 0;
			for (int i = 0; i <= 1500; i++) {
				int index = (i / num);
				String file = outputFile + "/" + "File_" + index;			
				if (!list.contains(file)) {
					list.add(file);
					map.put(file, idx);
					System.out.println(file);
					idx++;
				}								
			}
			
			FileWriter fwriter2 = new FileWriter(outputScript);
		    BufferedWriter out2 = new BufferedWriter(fwriter2);
		    
			Iterator itr = map.keySet().iterator();
			while (itr.hasNext()) {
				String file = (String)itr.next();
				out2.write(jumpjarPath + " -GeneratePubchemStructureDatabase " + file + " " + outputFolder + "\n");
			}
			out2.close();
			
			FileWriter[] fwriter = new FileWriter[idx];
		    BufferedWriter[] out = new BufferedWriter[idx];
		    for (int i = 0; i < idx; i++) {
		    	fwriter[i] = new FileWriter((String)list.get(i));
		    	out[i] = new BufferedWriter(fwriter[i]);		    	
		    }
		    
			//int num = new Integer(args[1]);
			//String[] files = new String[num];
			FileInputStream fstream = new FileInputStream(inputFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				String id = split[0];
				String formula = split[1];
				String inchikey = split[2];
				String inchi = split[3];
				String smile = split[4];
				String iupac = split[5];
				double calc_mass = ToolBox.getMonoisotopicMass(split[1]);
				int mass_category = new Integer((int)(calc_mass / num));
				if (mass_category > 1500 / num) {
					mass_category = 1500 / num;
				}
				String file = outputFile + "/" + "File_" + mass_category;
				idx = (Integer)map.get(file);
				out[idx].write(str + "\n");
				out[idx].flush();
				
				
			}
			in.close();
			
			
			for (int i = 0; i < idx; i++) {
		    	
		    	out[i].close();
		    	
		    }
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
