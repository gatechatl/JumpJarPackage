package METABOLOMIC_DATABASE.PUBCHEM;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;

public class CompareTargetDecoyNumber {

	
	public static void execute(String[] args) {
		
		try {
			
			HashMap target_map = new HashMap();
			HashMap decoy_map = new HashMap();
			
			String outputFile = args[2];
			FileWriter fwriter = new FileWriter(outputFile);
            BufferedWriter out = new BufferedWriter(fwriter);			
            
            String outputFile_filtered = args[3];
			FileWriter fwriter_filtered = new FileWriter(outputFile_filtered);
            BufferedWriter out_filtered = new BufferedWriter(fwriter_filtered);			            
            
			String fileName = args[0];
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				String target = split[0].split(":")[0];
				String decoy = split[1].split(":")[0];
				
				target_map.put(target, target);
				decoy_map.put(decoy, decoy);
			}
			in.close();
			
			HashMap overlap_map = new HashMap();
			int overlap = 0;
			Iterator itr = target_map.keySet().iterator();
			while (itr.hasNext()) {
				String target = (String)itr.next();
				if (decoy_map.containsKey(target)) {
					overlap++;
					overlap_map.put(target, target);
					//System.out.println(target);
					out.write(target + "\n");
				}
			}
			out.close();
			
			fileName = args[1];
			fstream = new FileInputStream(fileName);
			din = new DataInputStream(fstream);
			in = new BufferedReader(new InputStreamReader(din));
			
			int count = 0;
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				out_filtered.write(split[0]);
				for (int i = 1; i < split.length; i++) {
					String target = split[i].split("\\|")[0];
					String target_formula = target.split(":")[0];
					String decoy = split[i].split("\\|")[1];
					String decoy_formula = decoy.split(":")[0];
					if (!overlap_map.containsKey(target_formula) && !overlap_map.containsKey(decoy_formula)) {
						System.out.println("what:" + target_formula);
						out_filtered.write("\t" + split[i]);
					} else {
						count++;
					}
				}								
				out_filtered.write("\n");
			}
			in.close();
			out_filtered.close();
			System.out.println(count);
			System.out.println(target_map.size() + "\t" + overlap + "\t" + decoy_map.size());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
