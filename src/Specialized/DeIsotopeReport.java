package Specialized;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.ObjectInputStream.GetField;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Examine the JUMPf result and compare it
 * @author tshaw
 *
 */
public class DeIsotopeReport {

	public static void execute(String[] args) {
		
		try {
			HashMap jumpf_map_count = new HashMap();
			HashMap jumpf_peptide = new HashMap();
			HashMap jumpf_peptide_mono = new HashMap();
			HashMap jumpf_peptideinfo = new HashMap();
			HashMap jumpf_map = new HashMap();
			String fileName = args[0];
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (split.length > 9 && !split[0].equals("Peptides")) {
					
					String peptide = split[0];
					String scan = split[7];
					String mass = split[8];
					String charge = split[9];
					String monopeakloc = split[14];
					PEPTIDE_INFO pinfo = new PEPTIDE_INFO();
					pinfo.PEPTIDE = peptide;
					pinfo.SCAN = new Integer(scan);
					
					
					pinfo.CHARGE = new Integer(charge);
					pinfo.MASS = ((new Double(mass) * pinfo.CHARGE) - 1.0078250321) / pinfo.CHARGE + 1.0078250321;
					pinfo.MONOPEAKLOC = new Integer(monopeakloc);
					
					int tag = new Double(pinfo.MASS).intValue();
					
					jumpf_peptide.put(pinfo.PEPTIDE, 0);
					jumpf_peptide_mono.put(pinfo.PEPTIDE, 0);
					jumpf_peptideinfo.put(pinfo.PEPTIDE, str);
					
					jumpf_map.put(pinfo.SCAN + ":" + (tag - 3), pinfo);
					jumpf_map.put(pinfo.SCAN + ":" + (tag - 2), pinfo);
					jumpf_map.put(pinfo.SCAN + ":" + (tag - 1), pinfo);
					jumpf_map.put(pinfo.SCAN + ":" + tag, pinfo);
					jumpf_map.put(pinfo.SCAN + ":" + (tag + 1), pinfo);
					jumpf_map.put(pinfo.SCAN + ":" + (tag + 2), pinfo);
					jumpf_map.put(pinfo.SCAN + ":" + (tag + 3), pinfo);
					
					jumpf_map_count.put(pinfo.SCAN + ":" + tag, 0);
					
				}
			}
			in.close();
			
			HashMap deisotope = new HashMap();
			String search_mass = "";
			String fileName2 = args[1];
			FileInputStream fstream2 = new FileInputStream(fileName2);
			DataInputStream din2 = new DataInputStream(fstream2);
			BufferedReader in2 = new BufferedReader(new InputStreamReader(din2));
			while (in2.ready()) {
				String str = in2.readLine();
				if (str.contains("searching: ")) {
					search_mass = str.split("searching: ")[1].split("\t")[0];
				}
				if (str.contains("scan number: ")) {
					int scan_number = new Integer(str.split("scan number: ")[1]);
					str = in2.readLine();
					int peak = new Integer(str.split("match for peak: ")[1]);
					str = in2.readLine();
					int charge = new Integer(str.split("match for charge: ")[1]);
					str = in2.readLine();
					String[] split = str.split("\t");
					double score = new Double(split[0]);
					double mono_isotopic_mass = new Double(split[17]);
					ISOTOPE_HIT hit = new ISOTOPE_HIT();
					hit.CHARGE = charge;
					hit.MONO = peak;
					hit.SCAN = scan_number;
					hit.SCORE = score;
					hit.MASS = new Double(search_mass);
					hit.MONOISOTOPICMASS = mono_isotopic_mass;
					if (hit.CHARGE == 1) {
						hit.SCORE = hit.SCORE - 0.2;
					}
					if (deisotope.containsKey(scan_number + ":" + search_mass)) {
						ISOTOPE_HIT old_hit = (ISOTOPE_HIT)deisotope.get(scan_number + ":" + search_mass);
						if (hit.SCORE > old_hit.SCORE) {
							deisotope.put(scan_number + ":" + search_mass, hit);
						}
					} else {
						deisotope.put(scan_number + ":" + search_mass, hit);
					}
				}
			}
			in2.close();
			
			int score_90 = 0;
			int score_80 = 0;
			int score_70 = 0;
			int score_60 = 0;
			int score_50 = 0;
			
			int also_found_hit = 0;
			int also_found_hit_same_charge = 0;
			int also_found_hit_same_charge_mono = 0;
			Iterator itr = deisotope.keySet().iterator();
			while (itr.hasNext()) {
				String key = (String)itr.next();
				ISOTOPE_HIT hit = (ISOTOPE_HIT)deisotope.get(key);
				String jumpf_key = hit.SCAN + ":" + new Double(hit.MONOISOTOPICMASS).intValue();
				if (hit.SCORE >= 0.9) {
					score_90++;
				} else if (hit.SCORE >= 0.8) {
					score_80++;
				} else if (hit.SCORE >= 0.7) {
					score_70++;
				} else if (hit.SCORE >= 0.6) {
					score_60++;
				} else if (hit.SCORE >= 0.5) {
					score_50++;
				}
				if (hit.SCORE >= 0.3 && jumpf_map.containsKey(jumpf_key)) {
					PEPTIDE_INFO pinfo = (PEPTIDE_INFO)jumpf_map.get(jumpf_key);
					also_found_hit++;
					double diff = Math.abs(pinfo.MASS - hit.MONOISOTOPICMASS);
					if (pinfo.CHARGE == hit.CHARGE && diff < 0.005){
						System.out.println(pinfo.PEPTIDE + "\t" + pinfo.SCAN + "\t" + pinfo.CHARGE + "\t" + pinfo.MASS + "\t" + hit.MONOISOTOPICMASS);
						
						also_found_hit_same_charge++;
						
						int pepcount = (Integer)jumpf_peptide.get(pinfo.PEPTIDE);
						jumpf_peptide.put(pinfo.PEPTIDE, pepcount + 1);
					}
					if (pinfo.CHARGE == hit.CHARGE && diff < 0.005 && pinfo.MONOPEAKLOC == hit.MONO){
						System.out.println(pinfo.PEPTIDE + "\t" + pinfo.SCAN + "\t" + pinfo.CHARGE + "\t" + pinfo.MASS + "\t" + hit.MONOISOTOPICMASS);
						
						also_found_hit_same_charge_mono++;
						
						int pepcount = (Integer)jumpf_peptide_mono.get(pinfo.PEPTIDE);
						jumpf_peptide_mono.put(pinfo.PEPTIDE, pepcount + 1);
					}
				}							
			}
			
			String missedoutputFile = args[2];
			FileWriter fwriter = new FileWriter(missedoutputFile);
            BufferedWriter out = new BufferedWriter(fwriter);
			
            String missedoutputFile_mono = args[3];
			FileWriter fwriter2 = new FileWriter(missedoutputFile_mono);
            BufferedWriter out2 = new BufferedWriter(fwriter2);
			
			int countget = 0;
			Iterator itr2 = jumpf_peptide.keySet().iterator();
			while (itr2.hasNext()) {
				String pep = (String)itr2.next();
				if ((Integer)jumpf_peptide.get(pep) > 0) {
					countget++;
				} else {
					
					out.write(jumpf_peptideinfo.get(pep) + "\n");
				}
			}
			out.close();
			
			int countget_mono = 0;
			itr2 = jumpf_peptide_mono.keySet().iterator();
			while (itr2.hasNext()) {
				String pep = (String)itr2.next();
				if ((Integer)jumpf_peptide_mono.get(pep) > 0) {
					countget_mono++;
				} else {
					
					out2.write(jumpf_peptideinfo.get(pep) + "\n");
				}
			}
			out2.close();
			System.out.println("Total number of jumpf: " + (jumpf_map.size() / 7));
			System.out.println("Total number of hits: " + (deisotope.size()));
			System.out.println("Total number of hits above 0.9: " + score_90);
			System.out.println("Total number of hits above 0.8: " + score_80);
			System.out.println("Total number of hits above 0.7: " + score_70);
			System.out.println("Total number of hits above 0.6: " + score_60);
			System.out.println("Total number of hits above 0.5: " + score_50);
			
			System.out.println("Total number of hits in jumpf: " + also_found_hit);
			System.out.println("Total number of hits in jumpf same charge: " + also_found_hit_same_charge);
			System.out.println("Total number of hits in jumpf same charge uniq: " + countget);
			System.out.println("Total number of hits in jumpf same charge mono uniq: " + countget_mono);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	public static class ISOTOPE_HIT {
		double MASS = -1;		
		int SCAN = -1;
		int MONO = -1;
		int CHARGE = -1;
		double SCORE = -1;
		double MONOISOTOPICMASS = -1;
	}
	public static class PEPTIDE_INFO {
		String PEPTIDE = "";
		int SCAN = -1;
		double MASS = -1;
		int CHARGE = -1;
		int MONOPEAKLOC = -1;
	}
}
