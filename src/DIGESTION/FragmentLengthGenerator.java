package DIGESTION;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

public class FragmentLengthGenerator {

	public static String parameter_info() {
		return "[fragmentFile]";
	}
	public static void execute(String[] args) {
		
		try {		
			String fragmentFile = args[0];
			boolean fasta = false;
			HashMap frag_list = new HashMap();
			int index = 0;
			
			FileInputStream fstream = new FileInputStream(fragmentFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				if (str.contains(">")) {
					fasta = true;
				}
				if (fasta) {
					if (str.contains(">")) {
						index++;
					} else {
						if (frag_list.containsKey(index)) {
							String seq = (String)frag_list.get(index);
							seq += str;
							frag_list.put(index, seq);
						} else {
							frag_list.put(index, str);
						}
					}
				} else {					
					index++;
					frag_list.put(index, str);
				}
			}
			in.close();
			
			Iterator itr = frag_list.keySet().iterator();
			while (itr.hasNext()) {
				index = (Integer)itr.next();
				String seq = (String)frag_list.get(index);
				System.out.println(seq.length());
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
