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
import MISC.ToolBox;

public class CreateIndexFormulaOriginal {
	public static void execute(String[] args) {
		try {
			HashMap all_formula = new HashMap();
			HashMap[] map = new HashMap[151000];
			for (int i = 0; i < 151000; i++) {
				map[i] = new HashMap();
			}

			String inputFile = args[0];
			String folder = args[1];
			double min = new Double(args[2]);
            double max = new Double(args[3]);
            /*String otherReference = args[4];
            String referenceNames = args[5];
            String[] splitRef = otherReference.split(",");
            String[] splitRefName = referenceNames.split(",");
            HashMap[] maps = new HashMap[splitRef.length];
            for (int i = 0; i < maps.length; i++) {
            	maps[i] = new HashMap();
    			FileInputStream fstream = new FileInputStream(splitRef[i]);
    			DataInputStream din = new DataInputStream(fstream);
    			BufferedReader in = new BufferedReader(new InputStreamReader(din));
    			while (in.ready()) {
    				String str = in.readLine();
    				String name = ToolBox.standardize_name(str.trim());
    				maps[i].put(name, name);
    			}
    			in.close();
    			
            }
            */
			FileInputStream fstream = new FileInputStream(inputFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (split.length > 2) {
					for (int i = 2; i < split.length; i++) {
						String[] split2 = split[i].split(":");
						double mass = new Double(split2[1]);					
						int index = (int)(mass * 100);
						if (mass >= min && mass < max) {
							map[index].put(split[i], "");
						}
					}
				}
			}
			in.close();			
			
			File f = new File(folder);
			if (!f.exists()) {
				f.mkdir();
			}
			
			for (int i = 0; i < 151000; i++) {
				if (i >= min * 100 && i < max * 100) {
				FileWriter fwriter = new FileWriter(folder + "/Formula" + (i + 1) + ".txt"); 
				BufferedWriter out = new BufferedWriter(fwriter);
				
				Iterator itr = map[i].keySet().iterator();
				while (itr.hasNext()) {
					String str = (String)itr.next();
					out.write(str + "\n");
					/*String[] split = str.split(":");
					String formula = ToolBox.standardize_name(split[0]);
					String mass = split[1];
					String hits = "";
					for (int j = 0; j < maps.length; j++) {
						if (maps[j].containsKey(formula)) {
							hits += splitRefName[j] + ",";
						}
					}
					if (hits.equals("")) {
						out.write(formula + ":" + mass + "\n");
					} else {
						out.write(formula + ":" + mass + "\t" + hits + "\n");
					}*/
				}
				out.close();
				}
			}						
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

