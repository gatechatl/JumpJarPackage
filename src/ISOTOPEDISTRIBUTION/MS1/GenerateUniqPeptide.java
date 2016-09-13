package ISOTOPEDISTRIBUTION.MS1;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

/**
 * Small function for removing duplicated rows
 * @author tshaw
 */
public class GenerateUniqPeptide {

	public static void main(String[] args) {
		
		try {
			HashMap map = new HashMap();
			
			String fileName = args[0];
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (!map.containsKey(split[7])) {
					System.out.println(str);
				}
				map.put(split[7], split[7]);
			}
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
