package Specialized;

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

public class GeneratePubchemDecoyDatabase {

	public static void execute(String[] args) {
		try {
			HashMap all_formula = new HashMap();
			HashMap[] map = new HashMap[151000];
			for (int i = 0; i < 151000; i++) {
				map[i] = new HashMap();
			}
			//double min = new Double(args[2]);
            //double max = new Double(args[3]);
            String type = args[2];
            int type_index = 0;
            if (type.equals("NORM")) {
            	type_index = 0;
            } else if (type.equals("DECOY")) {
            	type_index = 1;
            }
            
			String inputFile = args[0];
			FileInputStream fstream = new FileInputStream(inputFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (split.length >= 2) {
					for (int i = 1; i < split.length; i++) {
						String[] split2 = split[i].split("\\|");
						
						
						String normal_formula = split2[0].split(":")[0];
						
						
						String text = split2[type_index];
						String[] split3 = text.split(":");
						String formula = split3[0];
						if (ToolBox.check_formula_valid_element(formula) && ToolBox.check_hydrogen_rule(normal_formula)) {
							double mass = new Double(split3[1]);
							
							int index = (int)(mass * 100);
							//if (mass > min && mass <= max) {
							if (index >= 151000) {
								index = 150999;
							}
							map[index].put(text, "");
							//}
						}
					}
				}
			}
			in.close();
			
			String folder = args[1];
			File f = new File(folder);
			if (!f.exists()) {
				f.mkdir();
			}
			
			for (int i = 0; i < 151000; i++) {
				//if (i > min * 100 && i <= max * 100) {
					FileWriter fwriter = new FileWriter(folder + "/Formula" + (i + 1) + ".txt"); 
					BufferedWriter out = new BufferedWriter(fwriter);
					
					Iterator itr = map[i].keySet().iterator();
					while (itr.hasNext()) {
						String str = (String)itr.next();
						out.write(str + "\n");
					}
					out.close();
				//}
			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
