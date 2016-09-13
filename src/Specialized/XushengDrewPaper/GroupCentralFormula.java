package Specialized.XushengDrewPaper;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * Identify the central equation
 * @author tshaw
 *
 */
public class GroupCentralFormula {

	public static void execute(String[] args) {
		
		try {
			
			String inputFile = args[0];
			String matrixFile = args[1];
			String origFile = args[2];
			// group each sample based on cluster put them into hashmap
			
			HashMap group = new HashMap();
			HashMap map = new HashMap();
			FileInputStream fstream = new FileInputStream(inputFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			in.readLine();
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split(" ");
				split[0] = split[0].replaceAll("\"", "");
				split[1] = split[1].replaceAll("\"", "");
				map.put(split[0], split[1]);
				group.put(split[1], split[1]);
			}
			in.close();
			
			HashMap smileID = new HashMap();
			int index = 1;
			fstream = new FileInputStream(origFile);
			din = new DataInputStream(fstream);
			in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				String id = "ID" + index;
				String smile = split[1];
				smile = smile.replaceAll("\\(", "(");
				smile = smile.replaceAll("\\)", ")");
				//System.out.println("perl JumpJar.pl -program=SmileFragmenter -Smile='" + smile + "' -mass=50 -depth=2 > fragment_result_" + id);
				smileID.put(id, smile);
				index++;
			}
			in.close();
			
			
			Iterator itr = group.keySet().iterator();
			while (itr.hasNext()) {								
				String current_sample_group = (String)itr.next();
				int total_features = 0;
				LinkedList all_sample = new LinkedList();
				// read the matrix table
				fstream = new FileInputStream(matrixFile);
				din = new DataInputStream(fstream);
				in = new BufferedReader(new InputStreamReader(din));
				in.readLine();
				while (in.ready()) {
					String str = in.readLine();
					String[] split = str.split("\t");
					total_features = split.length;
					String sample_group = (String)map.get(split[0]);
					if (sample_group.equals(current_sample_group)) {
						all_sample.add(str);
					}
				}
				in.close();
				
				
				HashMap count_col = new HashMap();
				Iterator itr4 = all_sample.iterator();
				while (itr4.hasNext()) {
					String str4 = (String)itr4.next();
					String[] split = str4.split("\t");
					for (int i = 1; i < total_features; i++) {
						if (split[i].equals("1")) {
							if (count_col.containsKey(i)) {
								int num = (Integer)count_col.get(i);
								count_col.put(i, (num + 1));				
							} else {
								count_col.put(i, 1);
							}
						}
					}
				}
				
				HashMap col_majority = new HashMap();
				for (int i = 1; i < total_features; i++) {
					if (count_col.containsKey(i)) {
						int num = (Integer)count_col.get(i);
						if (num * 2 >= all_sample.size()) {
							col_majority.put(i, new Double(num) / all_sample.size());
						}
					} else {
						//col_majority.put(i, 0);
					}
				}
				
				// compare each sample to each other for each sample
				HashMap score = new HashMap();
				double max_score = 0;
				String highest_sample = "";
				Iterator itr2 = all_sample.iterator();
				while (itr2.hasNext()) {
					String sample1 = (String)itr2.next();
					String[] split_sample1 = sample1.split("\t");
					double totalScore = 0;
					for (int i = 1; i < total_features; i++) {
						if (col_majority.containsKey(i) && split_sample1[i].equals("1")) {
							double value = (Double)col_majority.get(i);
							totalScore = totalScore + value;
							//totalScore++;
							
						} else if (col_majority.containsKey(i) && split_sample1[i].equals("0")) {
							double value = 1 - (Double)col_majority.get(i);
							totalScore = totalScore + value;
							//totalScore++;
							
						} else if (!col_majority.containsKey(i) && split_sample1[i].equals("0")){
							//double value = 1 - (Double)col_majority.get(i);
							//totalScore = totalScore + value;
							totalScore++;
						}
					}
					score.put(split_sample1[0], totalScore);
					if (max_score < totalScore) {
						if (!(((String)smileID.get(sample1.split("\t")[0])).contains("+") || ((String)smileID.get(sample1.split("\t")[0])).contains("."))) {
							max_score = totalScore;
							highest_sample = sample1.split("\t")[0];
						}
					} else if (max_score == totalScore && !((String)smileID.get(sample1.split("\t")[0])).contains("+")) {
						highest_sample = sample1.split("\t")[0];
					}
					//Iterator itr3 = all_sample.iterator();
					//while (itr3.hasNext()) {
						//String sample2 = (String)itr3.next();
						//String[] split_sample2 = sample2.split("\t");
						//compareString(sample1, sample2);
					//}					
				}
				System.out.println("Group: " + current_sample_group + " " + max_score + " " + highest_sample + " " + (String)smileID.get(highest_sample));								
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public static int compareString(String str1, String str2) {
		String[] split1 = str1.split("\t");
		String[] split2 = str2.split("\t");
		int count = 0;
		for (int i = 0; i < split1.length; i++) {
			if (split1[i].equals(split2[i])) {
				count++;
			}
		}
		return count;
	}
}
