package Specialized.XushengDrewPaper;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;

public class AppendData2FragmentCluster {
	public static void execute(String[] args) {
		try {
			
			
			String inputFile = args[0];
			String fragmentResult = args[1];
			String hmdbFile = args[2];
			String outputFile = args[3];
			
			FileWriter fwriter = new FileWriter(outputFile);
			BufferedWriter out = new BufferedWriter(fwriter);
			
			HashMap map = new HashMap();
			FileInputStream fstream = new FileInputStream(fragmentResult);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			in.readLine();
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split(" ");
				split[0] = split[0].replaceAll("\"", "");
				
				map.put(split[0], split[1]);
			}
			in.close();
			
			HashMap smileList = new HashMap();
			fstream = new FileInputStream(hmdbFile);
			din = new DataInputStream(fstream);
			in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (split.length > 4) {
					String smile = split[4];
					smile = smile.replaceAll("\\(", "\\\\(");
					smile = smile.replaceAll("\\)", "\\\\)");
					smileList.put(smile, smile);
				}
			}
			in.close();
			
			int index = 1;
			fstream = new FileInputStream(inputFile);
			din = new DataInputStream(fstream);
			in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				String smile = split[1];
				boolean hmdb_present = false;
				if (smileList.containsKey(smile)) {
					hmdb_present = true;
				}
				String id = "ID" + index;
				String cluster = (String)map.get(id);
				out.write(str + "\t" + cluster + "\t" + hmdb_present + "\n");
				index++;
			}
			in.close();
			out.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
