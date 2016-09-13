package MISC;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.ChemObject;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.formula.IsotopeContainer;
import org.openscience.cdk.formula.IsotopePattern;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import CDKFunction.IsotopePatternSimilarity;

public class ToolBox extends ChemObject {

	public static void main(String[] args) {
		try {
			System.out.println(check_formula_valid_element(standardize_name("C2H6OTl")));
			//convertInchi2SMILE("InChIKey=GPRLSGONYQIRFK-UHFFFAOYSA-N");
		} catch (Exception e) {
			e.printStackTrace();
		}
		/*String formula = "C38H26Cl4N11O6S4";
		System.out.println(check_hydrogen_rule(formula));
		System.exit(0);;
		System.out.println(standardize_name(formula));
		System.out.println(check_formula_valid_element(standardize_name("C9H11Br1Cl1S1P1NO2Cu1")));
		System.out.println(getMonoisotopicMass("C9H11NO2F1"));*/
	}
	
	public static String convertInchi2SMILE(String inchiStr) throws CDKException {
		InChIGeneratorFactory factory = new InChIGeneratorFactory();
		
		InChIToStructure structure = factory.getInChIToStructure(inchiStr, DefaultChemObjectBuilder.getInstance());
		IAtomContainer atom = structure.getAtomContainer();
		
		DefaultChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		SmilesParser sp = new SmilesParser(builder);
		
		// IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);
		// result = MolecularFormulaManipulator.getString(formula);
		SmilesGenerator generator = new SmilesGenerator();		
		generator.createSMILES((IMolecule)atom);
		//IsotopeContainer molecule = new IsotopeContainer();
		
		//generator.createSMILES()
		
		
		
		return inchiStr;		
		
	}
	
	public static boolean isNumeric(String str)
	{
	  return str.matches("-?\\d+(\\.\\d+)?");  //match a number with optional '-' and decimal.
	}
	/**
	 * Merge Resolution
	 * 
	 * @param mass1
	 *            other mass
	 * @param mass2
	 *            theoretical mass
	 * @param ppm_cutoff
	 * @return
	 */
	public static boolean within_distance(double mass1, double mass2,
			double ppm_cutoff) {
		
		double difference;
		if (mass1 >= mass2) {
			difference = mass1 - mass2;
		} else {
			difference = mass2 - mass1;
		}

		if (ppm_cutoff >= difference) {

			return true;
		}
		return false;
	}

	/**
	 * PPM = 1E6 * (Mass1 - Mass2) / Mass2
	 * 
	 * @param mass1
	 *            other mass
	 * @param mass2
	 *            theoretical mass
	 * @param ppm_cutoff
	 * @return
	 */
	public static boolean check_within_ppm(double mass1, double mass2,
			double ppm_cutoff) {
		double million = 1e6;
		double difference;
		if (mass1 >= mass2) {
			difference = mass1 - mass2;
		} else {
			difference = mass2 - mass1; 
		}
		double ppm = (million * difference) / ((mass2 + mass1) / 2);
		if (ppm_cutoff >= ppm) {
			return true;
		}
		return false;
	}
	public static double comparePeaks(LinkedList list1, LinkedList list2) {
		IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
		IsotopePattern pattern1 = getPeakInfo(list1);
		IsotopePattern pattern2 = getPeakInfo(list2);
		double similarity_score = similarity.compare(pattern1, pattern2);
		return similarity_score;
	}
	/**
	 * based on the linked list file generate isotopepattern object
	 * @param fileName
	 * @return
	 */
	public static IsotopePattern getPeakInfo(LinkedList list) {
		
		try {

			IsotopePattern result = new IsotopePattern();
			
			Iterator itr = list.iterator();
			while (itr.hasNext()) {
				String str = (String)itr.next();
				if (!str.trim().equals("")) {
					String[] split = str.split("\t");
					double mass = new Double(split[0]);
					double intensity = new Double(split[1]);
					result.addIsotope(new IsotopeContainer(mass, intensity));										
					//System.out.println`(mass + "\t" + intensity);
				}
			}
			
			return result;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public static String HillSystemOrder_Remove_H(String formula) {
		if (!isNumeric(formula.substring(formula.length() - 1, formula.length()))) {
			formula = formula + "1";			
		}
		int Ci_num = retrieve_num_element(formula, "Ci");
		int H_num = retrieve_num_element(formula, "H") - 1;
		int He_num = retrieve_num_element(formula, "He");
		int Li_num = retrieve_num_element(formula, "Li");
		int Be_num = retrieve_num_element(formula, "Be");
		int B_num = retrieve_num_element(formula, "B");
		int C_num = retrieve_num_element(formula, "C");
		int N_num = retrieve_num_element(formula, "N");
		int O_num = retrieve_num_element(formula, "O");
		int F_num = retrieve_num_element(formula, "F");
		int Ne_num = retrieve_num_element(formula, "Ne");
		int Na_num = retrieve_num_element(formula, "Na");
		int Mg_num = retrieve_num_element(formula, "Mg");
		int Al_num = retrieve_num_element(formula, "Al");
		int Si_num = retrieve_num_element(formula, "Si");
		int P_num = retrieve_num_element(formula, "P");
		int S_num = retrieve_num_element(formula, "S");
		int Cl_num = retrieve_num_element(formula, "Cl");
		int Ar_num = retrieve_num_element(formula, "Ar");
		int K_num = retrieve_num_element(formula, "K");
		int Ca_num = retrieve_num_element(formula, "Ca");
		int Sc_num = retrieve_num_element(formula, "Sc");
		int Ti_num = retrieve_num_element(formula, "Ti");
		int V_num = retrieve_num_element(formula, "V");
		int Cr_num = retrieve_num_element(formula, "Cr");
		int Mn_num = retrieve_num_element(formula, "Mn");
		int Fe_num = retrieve_num_element(formula, "Fe");
		int Co_num = retrieve_num_element(formula, "Co");
		int Ni_num = retrieve_num_element(formula, "Ni");
		int Cu_num = retrieve_num_element(formula, "Cu");
		int Zn_num = retrieve_num_element(formula, "Zn");
		int Ga_num = retrieve_num_element(formula, "Ga");
		int Ge_num = retrieve_num_element(formula, "Ge");
		int As_num = retrieve_num_element(formula, "As");
		int Se_num = retrieve_num_element(formula, "Se");
		int Br_num = retrieve_num_element(formula, "Br");
		int Kr_num = retrieve_num_element(formula, "Kr");
		int Rb_num = retrieve_num_element(formula, "Rb");
		int Sr_num = retrieve_num_element(formula, "Sr");
		int Y_num = retrieve_num_element(formula, "Y");
		int Zr_num = retrieve_num_element(formula, "Zr");
		int Nb_num = retrieve_num_element(formula, "Nb");
		int Mo_num = retrieve_num_element(formula, "Mo");
		int Tc_num = retrieve_num_element(formula, "Tc");
		int Ru_num = retrieve_num_element(formula, "Ru");
		int Rh_num = retrieve_num_element(formula, "Rh");
		int Pd_num = retrieve_num_element(formula, "Pd");
		int Ag_num = retrieve_num_element(formula, "Ag");
		int Cd_num = retrieve_num_element(formula, "Cd");
		int In_num = retrieve_num_element(formula, "In");
		int Sn_num = retrieve_num_element(formula, "Sn");
		int Sb_num = retrieve_num_element(formula, "Sb");
		int Te_num = retrieve_num_element(formula, "Te");
		int I_num = retrieve_num_element(formula, "I");
		int Xe_num = retrieve_num_element(formula, "Xe");
		int Cs_num = retrieve_num_element(formula, "Cs");
		int Ba_num = retrieve_num_element(formula, "Ba");
		int La_num = retrieve_num_element(formula, "La");
		int Ce_num = retrieve_num_element(formula, "Ce");
		int Pr_num = retrieve_num_element(formula, "Pr");
		int Nd_num = retrieve_num_element(formula, "Nd");
		int Pm_num = retrieve_num_element(formula, "Pm");
		int Sm_num = retrieve_num_element(formula, "Sm");
		int Eu_num = retrieve_num_element(formula, "Eu");
		int Gd_num = retrieve_num_element(formula, "Gd");
		int Tb_num = retrieve_num_element(formula, "Tb");
		int Dy_num = retrieve_num_element(formula, "Dy");
		int Ho_num = retrieve_num_element(formula, "Ho");
		int Er_num = retrieve_num_element(formula, "Er");
		int Tm_num = retrieve_num_element(formula, "Tm");
		int Yb_num = retrieve_num_element(formula, "Yb");
		int Lu_num = retrieve_num_element(formula, "Lu");
		int Hf_num = retrieve_num_element(formula, "Hf");
		int Ta_num = retrieve_num_element(formula, "Ta");
		int W_num = retrieve_num_element(formula, "W");
		int Re_num = retrieve_num_element(formula, "Re");
		int Os_num = retrieve_num_element(formula, "Os");
		int Ir_num = retrieve_num_element(formula, "Ir");
		int Pt_num = retrieve_num_element(formula, "Pt");
		int Au_num = retrieve_num_element(formula, "Au");
		int Hg_num = retrieve_num_element(formula, "Hg");
		int Tl_num = retrieve_num_element(formula, "Tl");
		int Pb_num = retrieve_num_element(formula, "Pb");
		int Bi_num = retrieve_num_element(formula, "Bi");
		int Po_num = retrieve_num_element(formula, "Po");
		int At_num = retrieve_num_element(formula, "At");
		int Rn_num = retrieve_num_element(formula, "Rn");
		int Fr_num = retrieve_num_element(formula, "Fr");
		int Ra_num = retrieve_num_element(formula, "Ra");
		int Ac_num = retrieve_num_element(formula, "Ac");
		int Th_num = retrieve_num_element(formula, "Th");
		int Pa_num = retrieve_num_element(formula, "Pa");
		int U_num = retrieve_num_element(formula, "U");
		int Np_num = retrieve_num_element(formula, "Np");
		int Pu_num = retrieve_num_element(formula, "Pu");
		int Am_num = retrieve_num_element(formula, "Am");
		int Cm_num = retrieve_num_element(formula, "Cm");
		int Bk_num = retrieve_num_element(formula, "Bk");
		int Cf_num = retrieve_num_element(formula, "Cf");
		int Es_num = retrieve_num_element(formula, "Es");
		int Fm_num = retrieve_num_element(formula, "Fm");
		int Md_num = retrieve_num_element(formula, "Md");
		int No_num = retrieve_num_element(formula, "No");
		int Lr_num = retrieve_num_element(formula, "Lr");
		int Rf_num = retrieve_num_element(formula, "Rf");
		int Db_num = retrieve_num_element(formula, "Db");
		int Sg_num = retrieve_num_element(formula, "Sg");
		int Bh_num = retrieve_num_element(formula, "Bh");
		int Hs_num = retrieve_num_element(formula, "Hs");
		int Mt_num = retrieve_num_element(formula, "Mt");
		int Ds_num = retrieve_num_element(formula, "Ds");
		int Rg_num = retrieve_num_element(formula, "Rg");
		int Cn_num = retrieve_num_element(formula, "Cn");
		int Uut_num = retrieve_num_element(formula, "Uut");
		int Fl_num = retrieve_num_element(formula, "Fl");
		int Uup_num = retrieve_num_element(formula, "Uup");
		int Lv_num = retrieve_num_element(formula, "Lv");
		int Uus_num = retrieve_num_element(formula, "Uus");
		int Uuo_num = retrieve_num_element(formula, "Uuo");
		String final_name = "";
		if (Ci_num > 0) {
			final_name += "Ci" + Ci_num;
		}
		if (C_num > 0) {
			final_name += "C" + C_num;
		}
		if (H_num > 0) {
			final_name += "H" + H_num;
		}
		if (Ar_num > 0) {
			final_name += "Ar" + Ar_num;
		}
		if (Al_num > 0) {
			final_name += "Al" + Al_num;
		}
		if (As_num > 0) {
			final_name += "As" + As_num;
		}
		if (Ag_num > 0) {
			final_name += "Ag" + Ag_num;
		}
		
		if (Br_num > 0) {
			final_name += "Br" + Br_num;
		}
		
		if (Be_num > 0) {
			final_name += "Be" + Be_num;
		}
		if (B_num > 0) {
			final_name += "B" + B_num;
		}
		if (Ca_num > 0) {
			final_name += "Ca" + Ca_num;
		}
		if (Cl_num > 0) {
			final_name += "Cl" + Cl_num;
		}
		if (Cr_num > 0) {
			final_name += "Cr" + Cr_num;
		}
		if (Co_num > 0) {
			final_name += "Co" + Co_num;
		}
		if (Cu_num > 0) {
			final_name += "Cu" + Cu_num;
		}
		
		if (F_num > 0) {
			final_name += "F" + F_num;
		}
		if (Fe_num > 0) {
			final_name += "Fe" + Fe_num;
		}
		if (Ga_num > 0) {
			final_name += "Ga" + Ga_num;
		}
		if (Ge_num > 0) {
			final_name += "Ge" + Ge_num;
		}
		
		if (He_num > 0) {
			final_name += "He" + He_num;
		}
		if (K_num > 0) {
			final_name += "K" + K_num;
		}
		if (Li_num > 0) {
			final_name += "Li" + Li_num;
		}
		if (Mg_num > 0) {
			final_name += "Mg" + Mg_num;
		}
		if (Mn_num > 0) {
			final_name += "Mn" + Mn_num;
		}
		if (Mo_num > 0) {
			final_name += "Mo" + Mo_num;
		}
		
		if (N_num > 0) {
			final_name += "N" + N_num;
		}
		if (Ni_num > 0) {
			final_name += "Ni" + Ni_num;
		}
		
		if (Ne_num > 0) {
			final_name += "Ne" + Ne_num;
		}
		if (Na_num > 0) {
			final_name += "Na" + Na_num;
		}
		if (Nb_num > 0) {
			final_name += "Nb" + Nb_num;
		}
		
		if (O_num > 0) {
			final_name += "O" + O_num;
		}		
		if (P_num > 0) {
			final_name += "P" + P_num;
		}						
		
		if (Si_num > 0) {
			final_name += "Si" + Si_num;
		}
		if (S_num > 0) {
			final_name += "S" + S_num;
		}
		
		
		if (Sc_num > 0) {
			final_name += "Sc" + Sc_num;
		}
		if (Se_num > 0) {
			final_name += "Se" + Se_num;
		}
		if (Sr_num > 0) {
			final_name += "Sr" + Sr_num;
		}
		
		if (Ti_num > 0) {
			final_name += "Ti" + Ti_num;
		}
		if (V_num > 0) {
			final_name += "V" + V_num;
		}
		if (Y_num > 0) {
			final_name += "Y" + Y_num;
		}
		if (Zr_num > 0) {
			final_name += "Zr" + Zr_num;
		}
		
		if (Zn_num > 0) {
			final_name += "Zn" + Zn_num;
		}
		
		if (Kr_num > 0) {
			final_name += "Kr" + Kr_num;
		}
		if (Rb_num > 0) {
			final_name += "Rb" + Rb_num;
		}
		if (Tc_num > 0) {
			final_name += "Tc" + Tc_num;
		}
		if (Ru_num > 0) {
			final_name += "Ru" + Ru_num;
		}
		if (Rh_num > 0) {
			final_name += "Rh" + Rh_num;
		}
		if (Pd_num > 0) {
			final_name += "Pd" + Pd_num;
		}
		if (Cd_num > 0) {
			final_name += "Cd" + Cd_num;
		}
		if (In_num > 0) {
			final_name += "In" + In_num;
		}
		if (Sn_num > 0) {
			final_name += "Sn" + Sn_num;
		}
		if (Sb_num > 0) {
			final_name += "Sb" + Sb_num;
		}
		if (Te_num > 0) {
			final_name += "Te" + Te_num;
		}
		if (I_num > 0) {
			final_name += "I" + I_num;
		}
		
		if (Cs_num > 0) {
			final_name += "Cs" + Cs_num;
		}
		if (Ba_num > 0) {
			final_name += "Ba" + Ba_num;
		}
		if (La_num > 0) {
			final_name += "La" + La_num;
		}
		if (Ce_num > 0) {
			final_name += "Ce" + Ce_num;
		}
		if (Pr_num > 0) {
			final_name += "Pr" + Pr_num;
		}
		if (Nd_num > 0) {
			final_name += "Nd" + Nd_num;
		}
		if (Pm_num > 0) {
			final_name += "Pm" + Pm_num;
		}
		if (Sm_num > 0) {
			final_name += "Sm" + Sm_num;
		}
		if (Eu_num > 0) {
			final_name += "Eu" + Eu_num;
		}
		if (Gd_num > 0) {
			final_name += "Gd" + Gd_num;
		}
		if (Tb_num > 0) {
			final_name += "Tb" + Tb_num;
		}
		if (Dy_num > 0) {
			final_name += "Dy" + Dy_num;
		}
		if (Ho_num > 0) {
			final_name += "Ho" + Ho_num;
		}
		if (Er_num > 0) {
			final_name += "Er" + Er_num;
		}
		if (Tm_num > 0) {
			final_name += "Tm" + Tm_num;
		}
		if (Yb_num > 0) {
			final_name += "Yb" + Yb_num;
		}
		if (Lu_num > 0) {
			final_name += "Lu" + Lu_num;
		}
		if (Hf_num > 0) {
			final_name += "Hf" + Hf_num;
		}
		if (Ta_num > 0) {
			final_name += "Ta" + Ta_num;
		}
		if (W_num > 0) {
			final_name += "W" + W_num;
		}
		if (Re_num > 0) {
			final_name += "Re" + Re_num;
		}
		if (Os_num > 0) {
			final_name += "Os" + Os_num;
		}
		if (Ir_num > 0) {
			final_name += "Ir" + Ir_num;
		}
		if (Pt_num > 0) {
			final_name += "Pt" + Pt_num;
		}
		if (Au_num > 0) {
			final_name += "Au" + Au_num;
		}
		if (Hg_num > 0) {
			final_name += "Hg" + Hg_num;
		}
		if (Tl_num > 0) {
			final_name += "Tl" + Tl_num;
		}
		if (Pb_num > 0) {
			final_name += "Pb" + Pb_num;
		}
		if (Bi_num > 0) {
			final_name += "Bi" + Bi_num;
		}
		if (Po_num > 0) {
			final_name += "Po" + Po_num;
		}
		if (At_num > 0) {
			final_name += "At" + At_num;
		}
		if (Rn_num > 0) {
			final_name += "Rn" + Rn_num;
		}
		if (Fr_num > 0) {
			final_name += "Fr" + Fr_num;
		}
		if (Ra_num > 0) {
			final_name += "Ra" + Ra_num;
		}
		if (Ac_num > 0) {
			final_name += "Ac" + Ac_num;
		}
		if (Th_num > 0) {
			final_name += "Th" + Th_num;
		}
		if (Pa_num > 0) {
			final_name += "Pa" + Pa_num;
		}
		if (U_num > 0) {
			final_name += "U" + U_num;
		}
		if (Np_num > 0) {
			final_name += "Np" + Np_num;
		}
		if (Pu_num > 0) {
			final_name += "Pu" + Pu_num;
		}
		if (Am_num > 0) {
			final_name += "Am" + Am_num;
		}
		if (Cm_num > 0) {
			final_name += "Cm" + Cm_num;
		}
		if (Bk_num > 0) {
			final_name += "Bk" + Bk_num;
		}
		if (Cf_num > 0) {
			final_name += "Cf" + Cf_num;
		}
		if (Es_num > 0) {
			final_name += "Es" + Es_num;
		}
		if (Fm_num > 0) {
			final_name += "Fm" + Fm_num;
		}
		if (Md_num > 0) {
			final_name += "Md" + Md_num;
		}
		if (No_num > 0) {
			final_name += "No" + No_num;
		}
		if (Lr_num > 0) {
			final_name += "Lr" + Lr_num;
		}
		if (Rf_num > 0) {
			final_name += "Rf" + Rf_num;
		}
		if (Db_num > 0) {
			final_name += "Db" + Db_num;
		}
		if (Sg_num > 0) {
			final_name += "Sg" + Sg_num;
		}
		if (Bh_num > 0) {
			final_name += "Bh" + Bh_num;
		}
		if (Hs_num > 0) {
			final_name += "Hs" + Hs_num;
		}
		if (Mt_num > 0) {
			final_name += "Mt" + Mt_num;
		}
		if (Ds_num > 0) {
			final_name += "Ds" + Ds_num;
		}
		if (Rg_num > 0) {
			final_name += "Rg" + Rg_num;
		}
		if (Cn_num > 0) {
			final_name += "Cn" + Cn_num;
		}
		if (Uut_num > 0) {
			final_name += "Uut" + Uut_num;
		}
		if (Fl_num > 0) {
			final_name += "Fl" + Fl_num;
		}
		if (Uup_num > 0) {
			final_name += "Uup" + Uup_num;
		}
		if (Lv_num > 0) {
			final_name += "Lv" + Lv_num;
		}
		if (Uus_num > 0) {
			final_name += "Uus" + Uus_num;
		}
		if (Uuo_num > 0) {
			final_name += "Uuo" + Uuo_num;
		}
		if (Xe_num > 0) {
			final_name += "Xe" + Xe_num;
		}
		return final_name;		
	}
	/**
	 * based on the input file, generate the IsotopePattern
	 * @param fileName
	 * @return
	 */
	public static IsotopePattern getPeakInfo(String fileName) {
		
		try {

			IsotopePattern result = new IsotopePattern();
			
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				if (!str.trim().equals("")) {
					String[] split = str.split("\t");
					double mass = new Double(split[0]);
					double intensity = new Double(split[1]);
					result.addIsotope(new IsotopeContainer(mass, intensity));										
					//System.out.println`(mass + "\t" + intensity);
				}
			}
			in.close();
			
			return result;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	public static String HillSystemOrder_ADD_H(String formula) {
		return HillSystemOrder_DECOY_H(formula);
	}
	public static String HillSystemOrder_DECOY_H(String formula) {
		if (!isNumeric(formula.substring(formula.length() - 1, formula.length()))) {
			formula = formula + "1";			
		}
		int Ci_num = retrieve_num_element(formula, "Ci");
		int H_num = retrieve_num_element(formula, "H") + 1;
		int He_num = retrieve_num_element(formula, "He");
		int Li_num = retrieve_num_element(formula, "Li");
		int Be_num = retrieve_num_element(formula, "Be");
		int B_num = retrieve_num_element(formula, "B");
		int C_num = retrieve_num_element(formula, "C");
		int N_num = retrieve_num_element(formula, "N");
		int O_num = retrieve_num_element(formula, "O");
		int F_num = retrieve_num_element(formula, "F");
		int Ne_num = retrieve_num_element(formula, "Ne");
		int Na_num = retrieve_num_element(formula, "Na");
		int Mg_num = retrieve_num_element(formula, "Mg");
		int Al_num = retrieve_num_element(formula, "Al");
		int Si_num = retrieve_num_element(formula, "Si");
		int P_num = retrieve_num_element(formula, "P");
		int S_num = retrieve_num_element(formula, "S");
		int Cl_num = retrieve_num_element(formula, "Cl");
		int Ar_num = retrieve_num_element(formula, "Ar");
		int K_num = retrieve_num_element(formula, "K");
		int Ca_num = retrieve_num_element(formula, "Ca");
		int Sc_num = retrieve_num_element(formula, "Sc");
		int Ti_num = retrieve_num_element(formula, "Ti");
		int V_num = retrieve_num_element(formula, "V");
		int Cr_num = retrieve_num_element(formula, "Cr");
		int Mn_num = retrieve_num_element(formula, "Mn");
		int Fe_num = retrieve_num_element(formula, "Fe");
		int Co_num = retrieve_num_element(formula, "Co");
		int Ni_num = retrieve_num_element(formula, "Ni");
		int Cu_num = retrieve_num_element(formula, "Cu");
		int Zn_num = retrieve_num_element(formula, "Zn");
		int Ga_num = retrieve_num_element(formula, "Ga");
		int Ge_num = retrieve_num_element(formula, "Ge");
		int As_num = retrieve_num_element(formula, "As");
		int Se_num = retrieve_num_element(formula, "Se");
		int Br_num = retrieve_num_element(formula, "Br");
		int Kr_num = retrieve_num_element(formula, "Kr");
		int Rb_num = retrieve_num_element(formula, "Rb");
		int Sr_num = retrieve_num_element(formula, "Sr");
		int Y_num = retrieve_num_element(formula, "Y");
		int Zr_num = retrieve_num_element(formula, "Zr");
		int Nb_num = retrieve_num_element(formula, "Nb");
		int Mo_num = retrieve_num_element(formula, "Mo");
		int Tc_num = retrieve_num_element(formula, "Tc");
		int Ru_num = retrieve_num_element(formula, "Ru");
		int Rh_num = retrieve_num_element(formula, "Rh");
		int Pd_num = retrieve_num_element(formula, "Pd");
		int Ag_num = retrieve_num_element(formula, "Ag");
		int Cd_num = retrieve_num_element(formula, "Cd");
		int In_num = retrieve_num_element(formula, "In");
		int Sn_num = retrieve_num_element(formula, "Sn");
		int Sb_num = retrieve_num_element(formula, "Sb");
		int Te_num = retrieve_num_element(formula, "Te");
		int I_num = retrieve_num_element(formula, "I");
		int Xe_num = retrieve_num_element(formula, "Xe");
		int Cs_num = retrieve_num_element(formula, "Cs");
		int Ba_num = retrieve_num_element(formula, "Ba");
		int La_num = retrieve_num_element(formula, "La");
		int Ce_num = retrieve_num_element(formula, "Ce");
		int Pr_num = retrieve_num_element(formula, "Pr");
		int Nd_num = retrieve_num_element(formula, "Nd");
		int Pm_num = retrieve_num_element(formula, "Pm");
		int Sm_num = retrieve_num_element(formula, "Sm");
		int Eu_num = retrieve_num_element(formula, "Eu");
		int Gd_num = retrieve_num_element(formula, "Gd");
		int Tb_num = retrieve_num_element(formula, "Tb");
		int Dy_num = retrieve_num_element(formula, "Dy");
		int Ho_num = retrieve_num_element(formula, "Ho");
		int Er_num = retrieve_num_element(formula, "Er");
		int Tm_num = retrieve_num_element(formula, "Tm");
		int Yb_num = retrieve_num_element(formula, "Yb");
		int Lu_num = retrieve_num_element(formula, "Lu");
		int Hf_num = retrieve_num_element(formula, "Hf");
		int Ta_num = retrieve_num_element(formula, "Ta");
		int W_num = retrieve_num_element(formula, "W");
		int Re_num = retrieve_num_element(formula, "Re");
		int Os_num = retrieve_num_element(formula, "Os");
		int Ir_num = retrieve_num_element(formula, "Ir");
		int Pt_num = retrieve_num_element(formula, "Pt");
		int Au_num = retrieve_num_element(formula, "Au");
		int Hg_num = retrieve_num_element(formula, "Hg");
		int Tl_num = retrieve_num_element(formula, "Tl");
		int Pb_num = retrieve_num_element(formula, "Pb");
		int Bi_num = retrieve_num_element(formula, "Bi");
		int Po_num = retrieve_num_element(formula, "Po");
		int At_num = retrieve_num_element(formula, "At");
		int Rn_num = retrieve_num_element(formula, "Rn");
		int Fr_num = retrieve_num_element(formula, "Fr");
		int Ra_num = retrieve_num_element(formula, "Ra");
		int Ac_num = retrieve_num_element(formula, "Ac");
		int Th_num = retrieve_num_element(formula, "Th");
		int Pa_num = retrieve_num_element(formula, "Pa");
		int U_num = retrieve_num_element(formula, "U");
		int Np_num = retrieve_num_element(formula, "Np");
		int Pu_num = retrieve_num_element(formula, "Pu");
		int Am_num = retrieve_num_element(formula, "Am");
		int Cm_num = retrieve_num_element(formula, "Cm");
		int Bk_num = retrieve_num_element(formula, "Bk");
		int Cf_num = retrieve_num_element(formula, "Cf");
		int Es_num = retrieve_num_element(formula, "Es");
		int Fm_num = retrieve_num_element(formula, "Fm");
		int Md_num = retrieve_num_element(formula, "Md");
		int No_num = retrieve_num_element(formula, "No");
		int Lr_num = retrieve_num_element(formula, "Lr");
		int Rf_num = retrieve_num_element(formula, "Rf");
		int Db_num = retrieve_num_element(formula, "Db");
		int Sg_num = retrieve_num_element(formula, "Sg");
		int Bh_num = retrieve_num_element(formula, "Bh");
		int Hs_num = retrieve_num_element(formula, "Hs");
		int Mt_num = retrieve_num_element(formula, "Mt");
		int Ds_num = retrieve_num_element(formula, "Ds");
		int Rg_num = retrieve_num_element(formula, "Rg");
		int Cn_num = retrieve_num_element(formula, "Cn");
		int Uut_num = retrieve_num_element(formula, "Uut");
		int Fl_num = retrieve_num_element(formula, "Fl");
		int Uup_num = retrieve_num_element(formula, "Uup");
		int Lv_num = retrieve_num_element(formula, "Lv");
		int Uus_num = retrieve_num_element(formula, "Uus");
		int Uuo_num = retrieve_num_element(formula, "Uuo");
		String final_name = "";
		if (Ci_num > 0) {
			final_name += "Ci" + Ci_num;
		}
		if (C_num > 0) {
			final_name += "C" + C_num;
		}
		if (H_num > 0) {
			final_name += "H" + H_num;
		}
		if (Ar_num > 0) {
			final_name += "Ar" + Ar_num;
		}
		if (Al_num > 0) {
			final_name += "Al" + Al_num;
		}
		if (As_num > 0) {
			final_name += "As" + As_num;
		}
		if (Ag_num > 0) {
			final_name += "Ag" + Ag_num;
		}
		
		if (Br_num > 0) {
			final_name += "Br" + Br_num;
		}
		
		if (Be_num > 0) {
			final_name += "Be" + Be_num;
		}
		if (B_num > 0) {
			final_name += "B" + B_num;
		}
		if (Ca_num > 0) {
			final_name += "Ca" + Ca_num;
		}
		if (Cl_num > 0) {
			final_name += "Cl" + Cl_num;
		}
		if (Cr_num > 0) {
			final_name += "Cr" + Cr_num;
		}
		if (Co_num > 0) {
			final_name += "Co" + Co_num;
		}
		if (Cu_num > 0) {
			final_name += "Cu" + Cu_num;
		}
		
		if (F_num > 0) {
			final_name += "F" + F_num;
		}
		if (Fe_num > 0) {
			final_name += "Fe" + Fe_num;
		}
		if (Ga_num > 0) {
			final_name += "Ga" + Ga_num;
		}
		if (Ge_num > 0) {
			final_name += "Ge" + Ge_num;
		}
		
		if (He_num > 0) {
			final_name += "He" + He_num;
		}
		if (K_num > 0) {
			final_name += "K" + K_num;
		}
		if (Li_num > 0) {
			final_name += "Li" + Li_num;
		}
		if (Mg_num > 0) {
			final_name += "Mg" + Mg_num;
		}
		if (Mn_num > 0) {
			final_name += "Mn" + Mn_num;
		}
		if (Mo_num > 0) {
			final_name += "Mo" + Mo_num;
		}
		
		if (N_num > 0) {
			final_name += "N" + N_num;
		}
		if (Ni_num > 0) {
			final_name += "Ni" + Ni_num;
		}
		
		if (Ne_num > 0) {
			final_name += "Ne" + Ne_num;
		}
		if (Na_num > 0) {
			final_name += "Na" + Na_num;
		}
		if (Nb_num > 0) {
			final_name += "Nb" + Nb_num;
		}
		
		if (O_num > 0) {
			final_name += "O" + O_num;
		}		
		if (P_num > 0) {
			final_name += "P" + P_num;
		}						
		
		if (Si_num > 0) {
			final_name += "Si" + Si_num;
		}
		if (S_num > 0) {
			final_name += "S" + S_num;
		}
		
		
		if (Sc_num > 0) {
			final_name += "Sc" + Sc_num;
		}
		if (Se_num > 0) {
			final_name += "Se" + Se_num;
		}
		if (Sr_num > 0) {
			final_name += "Sr" + Sr_num;
		}
		
		if (Ti_num > 0) {
			final_name += "Ti" + Ti_num;
		}
		if (V_num > 0) {
			final_name += "V" + V_num;
		}
		if (Y_num > 0) {
			final_name += "Y" + Y_num;
		}
		if (Zr_num > 0) {
			final_name += "Zr" + Zr_num;
		}
		
		if (Zn_num > 0) {
			final_name += "Zn" + Zn_num;
		}
		
		if (Kr_num > 0) {
			final_name += "Kr" + Kr_num;
		}
		if (Rb_num > 0) {
			final_name += "Rb" + Rb_num;
		}
		if (Tc_num > 0) {
			final_name += "Tc" + Tc_num;
		}
		if (Ru_num > 0) {
			final_name += "Ru" + Ru_num;
		}
		if (Rh_num > 0) {
			final_name += "Rh" + Rh_num;
		}
		if (Pd_num > 0) {
			final_name += "Pd" + Pd_num;
		}
		if (Cd_num > 0) {
			final_name += "Cd" + Cd_num;
		}
		if (In_num > 0) {
			final_name += "In" + In_num;
		}
		if (Sn_num > 0) {
			final_name += "Sn" + Sn_num;
		}
		if (Sb_num > 0) {
			final_name += "Sb" + Sb_num;
		}
		if (Te_num > 0) {
			final_name += "Te" + Te_num;
		}
		if (I_num > 0) {
			final_name += "I" + I_num;
		}
		
		if (Cs_num > 0) {
			final_name += "Cs" + Cs_num;
		}
		if (Ba_num > 0) {
			final_name += "Ba" + Ba_num;
		}
		if (La_num > 0) {
			final_name += "La" + La_num;
		}
		if (Ce_num > 0) {
			final_name += "Ce" + Ce_num;
		}
		if (Pr_num > 0) {
			final_name += "Pr" + Pr_num;
		}
		if (Nd_num > 0) {
			final_name += "Nd" + Nd_num;
		}
		if (Pm_num > 0) {
			final_name += "Pm" + Pm_num;
		}
		if (Sm_num > 0) {
			final_name += "Sm" + Sm_num;
		}
		if (Eu_num > 0) {
			final_name += "Eu" + Eu_num;
		}
		if (Gd_num > 0) {
			final_name += "Gd" + Gd_num;
		}
		if (Tb_num > 0) {
			final_name += "Tb" + Tb_num;
		}
		if (Dy_num > 0) {
			final_name += "Dy" + Dy_num;
		}
		if (Ho_num > 0) {
			final_name += "Ho" + Ho_num;
		}
		if (Er_num > 0) {
			final_name += "Er" + Er_num;
		}
		if (Tm_num > 0) {
			final_name += "Tm" + Tm_num;
		}
		if (Yb_num > 0) {
			final_name += "Yb" + Yb_num;
		}
		if (Lu_num > 0) {
			final_name += "Lu" + Lu_num;
		}
		if (Hf_num > 0) {
			final_name += "Hf" + Hf_num;
		}
		if (Ta_num > 0) {
			final_name += "Ta" + Ta_num;
		}
		if (W_num > 0) {
			final_name += "W" + W_num;
		}
		if (Re_num > 0) {
			final_name += "Re" + Re_num;
		}
		if (Os_num > 0) {
			final_name += "Os" + Os_num;
		}
		if (Ir_num > 0) {
			final_name += "Ir" + Ir_num;
		}
		if (Pt_num > 0) {
			final_name += "Pt" + Pt_num;
		}
		if (Au_num > 0) {
			final_name += "Au" + Au_num;
		}
		if (Hg_num > 0) {
			final_name += "Hg" + Hg_num;
		}
		if (Tl_num > 0) {
			final_name += "Tl" + Tl_num;
		}
		if (Pb_num > 0) {
			final_name += "Pb" + Pb_num;
		}
		if (Bi_num > 0) {
			final_name += "Bi" + Bi_num;
		}
		if (Po_num > 0) {
			final_name += "Po" + Po_num;
		}
		if (At_num > 0) {
			final_name += "At" + At_num;
		}
		if (Rn_num > 0) {
			final_name += "Rn" + Rn_num;
		}
		if (Fr_num > 0) {
			final_name += "Fr" + Fr_num;
		}
		if (Ra_num > 0) {
			final_name += "Ra" + Ra_num;
		}
		if (Ac_num > 0) {
			final_name += "Ac" + Ac_num;
		}
		if (Th_num > 0) {
			final_name += "Th" + Th_num;
		}
		if (Pa_num > 0) {
			final_name += "Pa" + Pa_num;
		}
		if (U_num > 0) {
			final_name += "U" + U_num;
		}
		if (Np_num > 0) {
			final_name += "Np" + Np_num;
		}
		if (Pu_num > 0) {
			final_name += "Pu" + Pu_num;
		}
		if (Am_num > 0) {
			final_name += "Am" + Am_num;
		}
		if (Cm_num > 0) {
			final_name += "Cm" + Cm_num;
		}
		if (Bk_num > 0) {
			final_name += "Bk" + Bk_num;
		}
		if (Cf_num > 0) {
			final_name += "Cf" + Cf_num;
		}
		if (Es_num > 0) {
			final_name += "Es" + Es_num;
		}
		if (Fm_num > 0) {
			final_name += "Fm" + Fm_num;
		}
		if (Md_num > 0) {
			final_name += "Md" + Md_num;
		}
		if (No_num > 0) {
			final_name += "No" + No_num;
		}
		if (Lr_num > 0) {
			final_name += "Lr" + Lr_num;
		}
		if (Rf_num > 0) {
			final_name += "Rf" + Rf_num;
		}
		if (Db_num > 0) {
			final_name += "Db" + Db_num;
		}
		if (Sg_num > 0) {
			final_name += "Sg" + Sg_num;
		}
		if (Bh_num > 0) {
			final_name += "Bh" + Bh_num;
		}
		if (Hs_num > 0) {
			final_name += "Hs" + Hs_num;
		}
		if (Mt_num > 0) {
			final_name += "Mt" + Mt_num;
		}
		if (Ds_num > 0) {
			final_name += "Ds" + Ds_num;
		}
		if (Rg_num > 0) {
			final_name += "Rg" + Rg_num;
		}
		if (Cn_num > 0) {
			final_name += "Cn" + Cn_num;
		}
		if (Uut_num > 0) {
			final_name += "Uut" + Uut_num;
		}
		if (Fl_num > 0) {
			final_name += "Fl" + Fl_num;
		}
		if (Uup_num > 0) {
			final_name += "Uup" + Uup_num;
		}
		if (Lv_num > 0) {
			final_name += "Lv" + Lv_num;
		}
		if (Uus_num > 0) {
			final_name += "Uus" + Uus_num;
		}
		if (Uuo_num > 0) {
			final_name += "Uuo" + Uuo_num;
		}
		if (Xe_num > 0) {
			final_name += "Xe" + Xe_num;
		}
		return final_name;		
	}
	public static String HillSystemOrder_DECOY_NH2(String formula) {
		if (!isNumeric(formula.substring(formula.length() - 1, formula.length()))) {
			formula = formula + "1";			
		}
		int Ci_num = retrieve_num_element(formula, "Ci");
		int H_num = retrieve_num_element(formula, "H") + 2;
		int He_num = retrieve_num_element(formula, "He");
		int Li_num = retrieve_num_element(formula, "Li");
		int Be_num = retrieve_num_element(formula, "Be");
		int B_num = retrieve_num_element(formula, "B");
		int C_num = retrieve_num_element(formula, "C");
		int N_num = retrieve_num_element(formula, "N") + 1;
		int O_num = retrieve_num_element(formula, "O");
		int F_num = retrieve_num_element(formula, "F");
		int Ne_num = retrieve_num_element(formula, "Ne");
		int Na_num = retrieve_num_element(formula, "Na");
		int Mg_num = retrieve_num_element(formula, "Mg");
		int Al_num = retrieve_num_element(formula, "Al");
		int Si_num = retrieve_num_element(formula, "Si");
		int P_num = retrieve_num_element(formula, "P");
		int S_num = retrieve_num_element(formula, "S");
		int Cl_num = retrieve_num_element(formula, "Cl");
		int Ar_num = retrieve_num_element(formula, "Ar");
		int K_num = retrieve_num_element(formula, "K");
		int Ca_num = retrieve_num_element(formula, "Ca");
		int Sc_num = retrieve_num_element(formula, "Sc");
		int Ti_num = retrieve_num_element(formula, "Ti");
		int V_num = retrieve_num_element(formula, "V");
		int Cr_num = retrieve_num_element(formula, "Cr");
		int Mn_num = retrieve_num_element(formula, "Mn");
		int Fe_num = retrieve_num_element(formula, "Fe");
		int Co_num = retrieve_num_element(formula, "Co");
		int Ni_num = retrieve_num_element(formula, "Ni");
		int Cu_num = retrieve_num_element(formula, "Cu");
		int Zn_num = retrieve_num_element(formula, "Zn");
		int Ga_num = retrieve_num_element(formula, "Ga");
		int Ge_num = retrieve_num_element(formula, "Ge");
		int As_num = retrieve_num_element(formula, "As");
		int Se_num = retrieve_num_element(formula, "Se");
		int Br_num = retrieve_num_element(formula, "Br");
		int Kr_num = retrieve_num_element(formula, "Kr");
		int Rb_num = retrieve_num_element(formula, "Rb");
		int Sr_num = retrieve_num_element(formula, "Sr");
		int Y_num = retrieve_num_element(formula, "Y");
		int Zr_num = retrieve_num_element(formula, "Zr");
		int Nb_num = retrieve_num_element(formula, "Nb");
		int Mo_num = retrieve_num_element(formula, "Mo");
		int Tc_num = retrieve_num_element(formula, "Tc");
		int Ru_num = retrieve_num_element(formula, "Ru");
		int Rh_num = retrieve_num_element(formula, "Rh");
		int Pd_num = retrieve_num_element(formula, "Pd");
		int Ag_num = retrieve_num_element(formula, "Ag");
		int Cd_num = retrieve_num_element(formula, "Cd");
		int In_num = retrieve_num_element(formula, "In");
		int Sn_num = retrieve_num_element(formula, "Sn");
		int Sb_num = retrieve_num_element(formula, "Sb");
		int Te_num = retrieve_num_element(formula, "Te");
		int I_num = retrieve_num_element(formula, "I");
		int Xe_num = retrieve_num_element(formula, "Xe");
		int Cs_num = retrieve_num_element(formula, "Cs");
		int Ba_num = retrieve_num_element(formula, "Ba");
		int La_num = retrieve_num_element(formula, "La");
		int Ce_num = retrieve_num_element(formula, "Ce");
		int Pr_num = retrieve_num_element(formula, "Pr");
		int Nd_num = retrieve_num_element(formula, "Nd");
		int Pm_num = retrieve_num_element(formula, "Pm");
		int Sm_num = retrieve_num_element(formula, "Sm");
		int Eu_num = retrieve_num_element(formula, "Eu");
		int Gd_num = retrieve_num_element(formula, "Gd");
		int Tb_num = retrieve_num_element(formula, "Tb");
		int Dy_num = retrieve_num_element(formula, "Dy");
		int Ho_num = retrieve_num_element(formula, "Ho");
		int Er_num = retrieve_num_element(formula, "Er");
		int Tm_num = retrieve_num_element(formula, "Tm");
		int Yb_num = retrieve_num_element(formula, "Yb");
		int Lu_num = retrieve_num_element(formula, "Lu");
		int Hf_num = retrieve_num_element(formula, "Hf");
		int Ta_num = retrieve_num_element(formula, "Ta");
		int W_num = retrieve_num_element(formula, "W");
		int Re_num = retrieve_num_element(formula, "Re");
		int Os_num = retrieve_num_element(formula, "Os");
		int Ir_num = retrieve_num_element(formula, "Ir");
		int Pt_num = retrieve_num_element(formula, "Pt");
		int Au_num = retrieve_num_element(formula, "Au");
		int Hg_num = retrieve_num_element(formula, "Hg");
		int Tl_num = retrieve_num_element(formula, "Tl");
		int Pb_num = retrieve_num_element(formula, "Pb");
		int Bi_num = retrieve_num_element(formula, "Bi");
		int Po_num = retrieve_num_element(formula, "Po");
		int At_num = retrieve_num_element(formula, "At");
		int Rn_num = retrieve_num_element(formula, "Rn");
		int Fr_num = retrieve_num_element(formula, "Fr");
		int Ra_num = retrieve_num_element(formula, "Ra");
		int Ac_num = retrieve_num_element(formula, "Ac");
		int Th_num = retrieve_num_element(formula, "Th");
		int Pa_num = retrieve_num_element(formula, "Pa");
		int U_num = retrieve_num_element(formula, "U");
		int Np_num = retrieve_num_element(formula, "Np");
		int Pu_num = retrieve_num_element(formula, "Pu");
		int Am_num = retrieve_num_element(formula, "Am");
		int Cm_num = retrieve_num_element(formula, "Cm");
		int Bk_num = retrieve_num_element(formula, "Bk");
		int Cf_num = retrieve_num_element(formula, "Cf");
		int Es_num = retrieve_num_element(formula, "Es");
		int Fm_num = retrieve_num_element(formula, "Fm");
		int Md_num = retrieve_num_element(formula, "Md");
		int No_num = retrieve_num_element(formula, "No");
		int Lr_num = retrieve_num_element(formula, "Lr");
		int Rf_num = retrieve_num_element(formula, "Rf");
		int Db_num = retrieve_num_element(formula, "Db");
		int Sg_num = retrieve_num_element(formula, "Sg");
		int Bh_num = retrieve_num_element(formula, "Bh");
		int Hs_num = retrieve_num_element(formula, "Hs");
		int Mt_num = retrieve_num_element(formula, "Mt");
		int Ds_num = retrieve_num_element(formula, "Ds");
		int Rg_num = retrieve_num_element(formula, "Rg");
		int Cn_num = retrieve_num_element(formula, "Cn");
		int Uut_num = retrieve_num_element(formula, "Uut");
		int Fl_num = retrieve_num_element(formula, "Fl");
		int Uup_num = retrieve_num_element(formula, "Uup");
		int Lv_num = retrieve_num_element(formula, "Lv");
		int Uus_num = retrieve_num_element(formula, "Uus");
		int Uuo_num = retrieve_num_element(formula, "Uuo");
		String final_name = "";
		if (Ci_num > 0) {
			final_name += "Ci" + Ci_num;
		}
		if (C_num > 0) {
			final_name += "C" + C_num;
		}
		if (H_num > 0) {
			final_name += "H" + H_num;
		}
		if (Ar_num > 0) {
			final_name += "Ar" + Ar_num;
		}
		if (Al_num > 0) {
			final_name += "Al" + Al_num;
		}
		if (As_num > 0) {
			final_name += "As" + As_num;
		}
		if (Ag_num > 0) {
			final_name += "Ag" + Ag_num;
		}
		
		if (Br_num > 0) {
			final_name += "Br" + Br_num;
		}
		
		if (Be_num > 0) {
			final_name += "Be" + Be_num;
		}
		if (B_num > 0) {
			final_name += "B" + B_num;
		}
		if (Ca_num > 0) {
			final_name += "Ca" + Ca_num;
		}
		if (Cl_num > 0) {
			final_name += "Cl" + Cl_num;
		}
		if (Cr_num > 0) {
			final_name += "Cr" + Cr_num;
		}
		if (Co_num > 0) {
			final_name += "Co" + Co_num;
		}
		if (Cu_num > 0) {
			final_name += "Cu" + Cu_num;
		}
		
		if (F_num > 0) {
			final_name += "F" + F_num;
		}
		if (Fe_num > 0) {
			final_name += "Fe" + Fe_num;
		}
		if (Ga_num > 0) {
			final_name += "Ga" + Ga_num;
		}
		if (Ge_num > 0) {
			final_name += "Ge" + Ge_num;
		}
		
		if (He_num > 0) {
			final_name += "He" + He_num;
		}
		if (K_num > 0) {
			final_name += "K" + K_num;
		}
		if (Li_num > 0) {
			final_name += "Li" + Li_num;
		}
		if (Mg_num > 0) {
			final_name += "Mg" + Mg_num;
		}
		if (Mn_num > 0) {
			final_name += "Mn" + Mn_num;
		}
		if (Mo_num > 0) {
			final_name += "Mo" + Mo_num;
		}
		
		if (N_num > 0) {
			final_name += "N" + N_num;
		}
		if (Ni_num > 0) {
			final_name += "Ni" + Ni_num;
		}
		
		if (Ne_num > 0) {
			final_name += "Ne" + Ne_num;
		}
		if (Na_num > 0) {
			final_name += "Na" + Na_num;
		}
		if (Nb_num > 0) {
			final_name += "Nb" + Nb_num;
		}
		
		if (O_num > 0) {
			final_name += "O" + O_num;
		}		
		if (P_num > 0) {
			final_name += "P" + P_num;
		}						
		
		if (Si_num > 0) {
			final_name += "Si" + Si_num;
		}
		if (S_num > 0) {
			final_name += "S" + S_num;
		}
		
		
		if (Sc_num > 0) {
			final_name += "Sc" + Sc_num;
		}
		if (Se_num > 0) {
			final_name += "Se" + Se_num;
		}
		if (Sr_num > 0) {
			final_name += "Sr" + Sr_num;
		}
		
		if (Ti_num > 0) {
			final_name += "Ti" + Ti_num;
		}
		if (V_num > 0) {
			final_name += "V" + V_num;
		}
		if (Y_num > 0) {
			final_name += "Y" + Y_num;
		}
		if (Zr_num > 0) {
			final_name += "Zr" + Zr_num;
		}
		
		if (Zn_num > 0) {
			final_name += "Zn" + Zn_num;
		}
		
		if (Kr_num > 0) {
			final_name += "Kr" + Kr_num;
		}
		if (Rb_num > 0) {
			final_name += "Rb" + Rb_num;
		}
		if (Tc_num > 0) {
			final_name += "Tc" + Tc_num;
		}
		if (Ru_num > 0) {
			final_name += "Ru" + Ru_num;
		}
		if (Rh_num > 0) {
			final_name += "Rh" + Rh_num;
		}
		if (Pd_num > 0) {
			final_name += "Pd" + Pd_num;
		}
		if (Cd_num > 0) {
			final_name += "Cd" + Cd_num;
		}
		if (In_num > 0) {
			final_name += "In" + In_num;
		}
		if (Sn_num > 0) {
			final_name += "Sn" + Sn_num;
		}
		if (Sb_num > 0) {
			final_name += "Sb" + Sb_num;
		}
		if (Te_num > 0) {
			final_name += "Te" + Te_num;
		}
		if (I_num > 0) {
			final_name += "I" + I_num;
		}
		
		if (Cs_num > 0) {
			final_name += "Cs" + Cs_num;
		}
		if (Ba_num > 0) {
			final_name += "Ba" + Ba_num;
		}
		if (La_num > 0) {
			final_name += "La" + La_num;
		}
		if (Ce_num > 0) {
			final_name += "Ce" + Ce_num;
		}
		if (Pr_num > 0) {
			final_name += "Pr" + Pr_num;
		}
		if (Nd_num > 0) {
			final_name += "Nd" + Nd_num;
		}
		if (Pm_num > 0) {
			final_name += "Pm" + Pm_num;
		}
		if (Sm_num > 0) {
			final_name += "Sm" + Sm_num;
		}
		if (Eu_num > 0) {
			final_name += "Eu" + Eu_num;
		}
		if (Gd_num > 0) {
			final_name += "Gd" + Gd_num;
		}
		if (Tb_num > 0) {
			final_name += "Tb" + Tb_num;
		}
		if (Dy_num > 0) {
			final_name += "Dy" + Dy_num;
		}
		if (Ho_num > 0) {
			final_name += "Ho" + Ho_num;
		}
		if (Er_num > 0) {
			final_name += "Er" + Er_num;
		}
		if (Tm_num > 0) {
			final_name += "Tm" + Tm_num;
		}
		if (Yb_num > 0) {
			final_name += "Yb" + Yb_num;
		}
		if (Lu_num > 0) {
			final_name += "Lu" + Lu_num;
		}
		if (Hf_num > 0) {
			final_name += "Hf" + Hf_num;
		}
		if (Ta_num > 0) {
			final_name += "Ta" + Ta_num;
		}
		if (W_num > 0) {
			final_name += "W" + W_num;
		}
		if (Re_num > 0) {
			final_name += "Re" + Re_num;
		}
		if (Os_num > 0) {
			final_name += "Os" + Os_num;
		}
		if (Ir_num > 0) {
			final_name += "Ir" + Ir_num;
		}
		if (Pt_num > 0) {
			final_name += "Pt" + Pt_num;
		}
		if (Au_num > 0) {
			final_name += "Au" + Au_num;
		}
		if (Hg_num > 0) {
			final_name += "Hg" + Hg_num;
		}
		if (Tl_num > 0) {
			final_name += "Tl" + Tl_num;
		}
		if (Pb_num > 0) {
			final_name += "Pb" + Pb_num;
		}
		if (Bi_num > 0) {
			final_name += "Bi" + Bi_num;
		}
		if (Po_num > 0) {
			final_name += "Po" + Po_num;
		}
		if (At_num > 0) {
			final_name += "At" + At_num;
		}
		if (Rn_num > 0) {
			final_name += "Rn" + Rn_num;
		}
		if (Fr_num > 0) {
			final_name += "Fr" + Fr_num;
		}
		if (Ra_num > 0) {
			final_name += "Ra" + Ra_num;
		}
		if (Ac_num > 0) {
			final_name += "Ac" + Ac_num;
		}
		if (Th_num > 0) {
			final_name += "Th" + Th_num;
		}
		if (Pa_num > 0) {
			final_name += "Pa" + Pa_num;
		}
		if (U_num > 0) {
			final_name += "U" + U_num;
		}
		if (Np_num > 0) {
			final_name += "Np" + Np_num;
		}
		if (Pu_num > 0) {
			final_name += "Pu" + Pu_num;
		}
		if (Am_num > 0) {
			final_name += "Am" + Am_num;
		}
		if (Cm_num > 0) {
			final_name += "Cm" + Cm_num;
		}
		if (Bk_num > 0) {
			final_name += "Bk" + Bk_num;
		}
		if (Cf_num > 0) {
			final_name += "Cf" + Cf_num;
		}
		if (Es_num > 0) {
			final_name += "Es" + Es_num;
		}
		if (Fm_num > 0) {
			final_name += "Fm" + Fm_num;
		}
		if (Md_num > 0) {
			final_name += "Md" + Md_num;
		}
		if (No_num > 0) {
			final_name += "No" + No_num;
		}
		if (Lr_num > 0) {
			final_name += "Lr" + Lr_num;
		}
		if (Rf_num > 0) {
			final_name += "Rf" + Rf_num;
		}
		if (Db_num > 0) {
			final_name += "Db" + Db_num;
		}
		if (Sg_num > 0) {
			final_name += "Sg" + Sg_num;
		}
		if (Bh_num > 0) {
			final_name += "Bh" + Bh_num;
		}
		if (Hs_num > 0) {
			final_name += "Hs" + Hs_num;
		}
		if (Mt_num > 0) {
			final_name += "Mt" + Mt_num;
		}
		if (Ds_num > 0) {
			final_name += "Ds" + Ds_num;
		}
		if (Rg_num > 0) {
			final_name += "Rg" + Rg_num;
		}
		if (Cn_num > 0) {
			final_name += "Cn" + Cn_num;
		}
		if (Uut_num > 0) {
			final_name += "Uut" + Uut_num;
		}
		if (Fl_num > 0) {
			final_name += "Fl" + Fl_num;
		}
		if (Uup_num > 0) {
			final_name += "Uup" + Uup_num;
		}
		if (Lv_num > 0) {
			final_name += "Lv" + Lv_num;
		}
		if (Uus_num > 0) {
			final_name += "Uus" + Uus_num;
		}
		if (Uuo_num > 0) {
			final_name += "Uuo" + Uuo_num;
		}
		if (Xe_num > 0) {
			final_name += "Xe" + Xe_num;
		}
		return final_name;		
	}
	public static String HillSystemOrder(String formula) {
		if (!isNumeric(formula.substring(formula.length() - 1, formula.length()))) {
			formula = formula + "1";			
		}
		int Ci_num = retrieve_num_element(formula, "Ci");
		int H_num = retrieve_num_element(formula, "H");
		int He_num = retrieve_num_element(formula, "He");
		int Li_num = retrieve_num_element(formula, "Li");
		int Be_num = retrieve_num_element(formula, "Be");
		int B_num = retrieve_num_element(formula, "B");
		int C_num = retrieve_num_element(formula, "C");
		int N_num = retrieve_num_element(formula, "N");
		int O_num = retrieve_num_element(formula, "O");
		int F_num = retrieve_num_element(formula, "F");
		int Ne_num = retrieve_num_element(formula, "Ne");
		int Na_num = retrieve_num_element(formula, "Na");
		int Mg_num = retrieve_num_element(formula, "Mg");
		int Al_num = retrieve_num_element(formula, "Al");
		int Si_num = retrieve_num_element(formula, "Si");
		int P_num = retrieve_num_element(formula, "P");
		int S_num = retrieve_num_element(formula, "S");
		int Cl_num = retrieve_num_element(formula, "Cl");
		int Ar_num = retrieve_num_element(formula, "Ar");
		int K_num = retrieve_num_element(formula, "K");
		int Ca_num = retrieve_num_element(formula, "Ca");
		int Sc_num = retrieve_num_element(formula, "Sc");
		int Ti_num = retrieve_num_element(formula, "Ti");
		int V_num = retrieve_num_element(formula, "V");
		int Cr_num = retrieve_num_element(formula, "Cr");
		int Mn_num = retrieve_num_element(formula, "Mn");
		int Fe_num = retrieve_num_element(formula, "Fe");
		int Co_num = retrieve_num_element(formula, "Co");
		int Ni_num = retrieve_num_element(formula, "Ni");
		int Cu_num = retrieve_num_element(formula, "Cu");
		int Zn_num = retrieve_num_element(formula, "Zn");
		int Ga_num = retrieve_num_element(formula, "Ga");
		int Ge_num = retrieve_num_element(formula, "Ge");
		int As_num = retrieve_num_element(formula, "As");
		int Se_num = retrieve_num_element(formula, "Se");
		int Br_num = retrieve_num_element(formula, "Br");
		int Kr_num = retrieve_num_element(formula, "Kr");
		int Rb_num = retrieve_num_element(formula, "Rb");
		int Sr_num = retrieve_num_element(formula, "Sr");
		int Y_num = retrieve_num_element(formula, "Y");
		int Zr_num = retrieve_num_element(formula, "Zr");
		int Nb_num = retrieve_num_element(formula, "Nb");
		int Mo_num = retrieve_num_element(formula, "Mo");
		int Tc_num = retrieve_num_element(formula, "Tc");
		int Ru_num = retrieve_num_element(formula, "Ru");
		int Rh_num = retrieve_num_element(formula, "Rh");
		int Pd_num = retrieve_num_element(formula, "Pd");
		int Ag_num = retrieve_num_element(formula, "Ag");
		int Cd_num = retrieve_num_element(formula, "Cd");
		int In_num = retrieve_num_element(formula, "In");
		int Sn_num = retrieve_num_element(formula, "Sn");
		int Sb_num = retrieve_num_element(formula, "Sb");
		int Te_num = retrieve_num_element(formula, "Te");
		int I_num = retrieve_num_element(formula, "I");
		int Xe_num = retrieve_num_element(formula, "Xe");
		int Cs_num = retrieve_num_element(formula, "Cs");
		int Ba_num = retrieve_num_element(formula, "Ba");
		int La_num = retrieve_num_element(formula, "La");
		int Ce_num = retrieve_num_element(formula, "Ce");
		int Pr_num = retrieve_num_element(formula, "Pr");
		int Nd_num = retrieve_num_element(formula, "Nd");
		int Pm_num = retrieve_num_element(formula, "Pm");
		int Sm_num = retrieve_num_element(formula, "Sm");
		int Eu_num = retrieve_num_element(formula, "Eu");
		int Gd_num = retrieve_num_element(formula, "Gd");
		int Tb_num = retrieve_num_element(formula, "Tb");
		int Dy_num = retrieve_num_element(formula, "Dy");
		int Ho_num = retrieve_num_element(formula, "Ho");
		int Er_num = retrieve_num_element(formula, "Er");
		int Tm_num = retrieve_num_element(formula, "Tm");
		int Yb_num = retrieve_num_element(formula, "Yb");
		int Lu_num = retrieve_num_element(formula, "Lu");
		int Hf_num = retrieve_num_element(formula, "Hf");
		int Ta_num = retrieve_num_element(formula, "Ta");
		int W_num = retrieve_num_element(formula, "W");
		int Re_num = retrieve_num_element(formula, "Re");
		int Os_num = retrieve_num_element(formula, "Os");
		int Ir_num = retrieve_num_element(formula, "Ir");
		int Pt_num = retrieve_num_element(formula, "Pt");
		int Au_num = retrieve_num_element(formula, "Au");
		int Hg_num = retrieve_num_element(formula, "Hg");
		int Tl_num = retrieve_num_element(formula, "Tl");
		int Pb_num = retrieve_num_element(formula, "Pb");
		int Bi_num = retrieve_num_element(formula, "Bi");
		int Po_num = retrieve_num_element(formula, "Po");
		int At_num = retrieve_num_element(formula, "At");
		int Rn_num = retrieve_num_element(formula, "Rn");
		int Fr_num = retrieve_num_element(formula, "Fr");
		int Ra_num = retrieve_num_element(formula, "Ra");
		int Ac_num = retrieve_num_element(formula, "Ac");
		int Th_num = retrieve_num_element(formula, "Th");
		int Pa_num = retrieve_num_element(formula, "Pa");
		int U_num = retrieve_num_element(formula, "U");
		int Np_num = retrieve_num_element(formula, "Np");
		int Pu_num = retrieve_num_element(formula, "Pu");
		int Am_num = retrieve_num_element(formula, "Am");
		int Cm_num = retrieve_num_element(formula, "Cm");
		int Bk_num = retrieve_num_element(formula, "Bk");
		int Cf_num = retrieve_num_element(formula, "Cf");
		int Es_num = retrieve_num_element(formula, "Es");
		int Fm_num = retrieve_num_element(formula, "Fm");
		int Md_num = retrieve_num_element(formula, "Md");
		int No_num = retrieve_num_element(formula, "No");
		int Lr_num = retrieve_num_element(formula, "Lr");
		int Rf_num = retrieve_num_element(formula, "Rf");
		int Db_num = retrieve_num_element(formula, "Db");
		int Sg_num = retrieve_num_element(formula, "Sg");
		int Bh_num = retrieve_num_element(formula, "Bh");
		int Hs_num = retrieve_num_element(formula, "Hs");
		int Mt_num = retrieve_num_element(formula, "Mt");
		int Ds_num = retrieve_num_element(formula, "Ds");
		int Rg_num = retrieve_num_element(formula, "Rg");
		int Cn_num = retrieve_num_element(formula, "Cn");
		int Uut_num = retrieve_num_element(formula, "Uut");
		int Fl_num = retrieve_num_element(formula, "Fl");
		int Uup_num = retrieve_num_element(formula, "Uup");
		int Lv_num = retrieve_num_element(formula, "Lv");
		int Uus_num = retrieve_num_element(formula, "Uus");
		int Uuo_num = retrieve_num_element(formula, "Uuo");
		String final_name = "";
		if (Ci_num > 0) {
			final_name += "Ci" + Ci_num;
		}
		if (C_num > 0) {
			final_name += "C" + C_num;
		}
		if (H_num > 0) {
			final_name += "H" + H_num;
		}
		if (Ar_num > 0) {
			final_name += "Ar" + Ar_num;
		}
		if (Al_num > 0) {
			final_name += "Al" + Al_num;
		}
		if (As_num > 0) {
			final_name += "As" + As_num;
		}
		if (Ag_num > 0) {
			final_name += "Ag" + Ag_num;
		}
		
		if (Br_num > 0) {
			final_name += "Br" + Br_num;
		}
		
		if (Be_num > 0) {
			final_name += "Be" + Be_num;
		}
		if (B_num > 0) {
			final_name += "B" + B_num;
		}
		if (Ca_num > 0) {
			final_name += "Ca" + Ca_num;
		}
		if (Cl_num > 0) {
			final_name += "Cl" + Cl_num;
		}
		if (Cr_num > 0) {
			final_name += "Cr" + Cr_num;
		}
		if (Co_num > 0) {
			final_name += "Co" + Co_num;
		}
		if (Cu_num > 0) {
			final_name += "Cu" + Cu_num;
		}
		
		if (F_num > 0) {
			final_name += "F" + F_num;
		}
		if (Fe_num > 0) {
			final_name += "Fe" + Fe_num;
		}
		if (Ga_num > 0) {
			final_name += "Ga" + Ga_num;
		}
		if (Ge_num > 0) {
			final_name += "Ge" + Ge_num;
		}
		
		if (He_num > 0) {
			final_name += "He" + He_num;
		}
		if (K_num > 0) {
			final_name += "K" + K_num;
		}
		if (Li_num > 0) {
			final_name += "Li" + Li_num;
		}
		if (Mg_num > 0) {
			final_name += "Mg" + Mg_num;
		}
		if (Mn_num > 0) {
			final_name += "Mn" + Mn_num;
		}
		if (Mo_num > 0) {
			final_name += "Mo" + Mo_num;
		}
		
		if (N_num > 0) {
			final_name += "N" + N_num;
		}
		if (Ni_num > 0) {
			final_name += "Ni" + Ni_num;
		}
		
		if (Ne_num > 0) {
			final_name += "Ne" + Ne_num;
		}
		if (Na_num > 0) {
			final_name += "Na" + Na_num;
		}
		if (Nb_num > 0) {
			final_name += "Nb" + Nb_num;
		}
		
		if (O_num > 0) {
			final_name += "O" + O_num;
		}		
		if (P_num > 0) {
			final_name += "P" + P_num;
		}						
		
		if (Si_num > 0) {
			final_name += "Si" + Si_num;
		}
		if (S_num > 0) {
			final_name += "S" + S_num;
		}
		
		
		if (Sc_num > 0) {
			final_name += "Sc" + Sc_num;
		}
		if (Se_num > 0) {
			final_name += "Se" + Se_num;
		}
		if (Sr_num > 0) {
			final_name += "Sr" + Sr_num;
		}
		
		if (Ti_num > 0) {
			final_name += "Ti" + Ti_num;
		}
		if (V_num > 0) {
			final_name += "V" + V_num;
		}
		if (Y_num > 0) {
			final_name += "Y" + Y_num;
		}
		if (Zr_num > 0) {
			final_name += "Zr" + Zr_num;
		}
		
		if (Zn_num > 0) {
			final_name += "Zn" + Zn_num;
		}
		
		if (Kr_num > 0) {
			final_name += "Kr" + Kr_num;
		}
		if (Rb_num > 0) {
			final_name += "Rb" + Rb_num;
		}
		if (Tc_num > 0) {
			final_name += "Tc" + Tc_num;
		}
		if (Ru_num > 0) {
			final_name += "Ru" + Ru_num;
		}
		if (Rh_num > 0) {
			final_name += "Rh" + Rh_num;
		}
		if (Pd_num > 0) {
			final_name += "Pd" + Pd_num;
		}
		if (Cd_num > 0) {
			final_name += "Cd" + Cd_num;
		}
		if (In_num > 0) {
			final_name += "In" + In_num;
		}
		if (Sn_num > 0) {
			final_name += "Sn" + Sn_num;
		}
		if (Sb_num > 0) {
			final_name += "Sb" + Sb_num;
		}
		if (Te_num > 0) {
			final_name += "Te" + Te_num;
		}
		if (I_num > 0) {
			final_name += "I" + I_num;
		}
		
		if (Cs_num > 0) {
			final_name += "Cs" + Cs_num;
		}
		if (Ba_num > 0) {
			final_name += "Ba" + Ba_num;
		}
		if (La_num > 0) {
			final_name += "La" + La_num;
		}
		if (Ce_num > 0) {
			final_name += "Ce" + Ce_num;
		}
		if (Pr_num > 0) {
			final_name += "Pr" + Pr_num;
		}
		if (Nd_num > 0) {
			final_name += "Nd" + Nd_num;
		}
		if (Pm_num > 0) {
			final_name += "Pm" + Pm_num;
		}
		if (Sm_num > 0) {
			final_name += "Sm" + Sm_num;
		}
		if (Eu_num > 0) {
			final_name += "Eu" + Eu_num;
		}
		if (Gd_num > 0) {
			final_name += "Gd" + Gd_num;
		}
		if (Tb_num > 0) {
			final_name += "Tb" + Tb_num;
		}
		if (Dy_num > 0) {
			final_name += "Dy" + Dy_num;
		}
		if (Ho_num > 0) {
			final_name += "Ho" + Ho_num;
		}
		if (Er_num > 0) {
			final_name += "Er" + Er_num;
		}
		if (Tm_num > 0) {
			final_name += "Tm" + Tm_num;
		}
		if (Yb_num > 0) {
			final_name += "Yb" + Yb_num;
		}
		if (Lu_num > 0) {
			final_name += "Lu" + Lu_num;
		}
		if (Hf_num > 0) {
			final_name += "Hf" + Hf_num;
		}
		if (Ta_num > 0) {
			final_name += "Ta" + Ta_num;
		}
		if (W_num > 0) {
			final_name += "W" + W_num;
		}
		if (Re_num > 0) {
			final_name += "Re" + Re_num;
		}
		if (Os_num > 0) {
			final_name += "Os" + Os_num;
		}
		if (Ir_num > 0) {
			final_name += "Ir" + Ir_num;
		}
		if (Pt_num > 0) {
			final_name += "Pt" + Pt_num;
		}
		if (Au_num > 0) {
			final_name += "Au" + Au_num;
		}
		if (Hg_num > 0) {
			final_name += "Hg" + Hg_num;
		}
		if (Tl_num > 0) {
			final_name += "Tl" + Tl_num;
		}
		if (Pb_num > 0) {
			final_name += "Pb" + Pb_num;
		}
		if (Bi_num > 0) {
			final_name += "Bi" + Bi_num;
		}
		if (Po_num > 0) {
			final_name += "Po" + Po_num;
		}
		if (At_num > 0) {
			final_name += "At" + At_num;
		}
		if (Rn_num > 0) {
			final_name += "Rn" + Rn_num;
		}
		if (Fr_num > 0) {
			final_name += "Fr" + Fr_num;
		}
		if (Ra_num > 0) {
			final_name += "Ra" + Ra_num;
		}
		if (Ac_num > 0) {
			final_name += "Ac" + Ac_num;
		}
		if (Th_num > 0) {
			final_name += "Th" + Th_num;
		}
		if (Pa_num > 0) {
			final_name += "Pa" + Pa_num;
		}
		if (U_num > 0) {
			final_name += "U" + U_num;
		}
		if (Np_num > 0) {
			final_name += "Np" + Np_num;
		}
		if (Pu_num > 0) {
			final_name += "Pu" + Pu_num;
		}
		if (Am_num > 0) {
			final_name += "Am" + Am_num;
		}
		if (Cm_num > 0) {
			final_name += "Cm" + Cm_num;
		}
		if (Bk_num > 0) {
			final_name += "Bk" + Bk_num;
		}
		if (Cf_num > 0) {
			final_name += "Cf" + Cf_num;
		}
		if (Es_num > 0) {
			final_name += "Es" + Es_num;
		}
		if (Fm_num > 0) {
			final_name += "Fm" + Fm_num;
		}
		if (Md_num > 0) {
			final_name += "Md" + Md_num;
		}
		if (No_num > 0) {
			final_name += "No" + No_num;
		}
		if (Lr_num > 0) {
			final_name += "Lr" + Lr_num;
		}
		if (Rf_num > 0) {
			final_name += "Rf" + Rf_num;
		}
		if (Db_num > 0) {
			final_name += "Db" + Db_num;
		}
		if (Sg_num > 0) {
			final_name += "Sg" + Sg_num;
		}
		if (Bh_num > 0) {
			final_name += "Bh" + Bh_num;
		}
		if (Hs_num > 0) {
			final_name += "Hs" + Hs_num;
		}
		if (Mt_num > 0) {
			final_name += "Mt" + Mt_num;
		}
		if (Ds_num > 0) {
			final_name += "Ds" + Ds_num;
		}
		if (Rg_num > 0) {
			final_name += "Rg" + Rg_num;
		}
		if (Cn_num > 0) {
			final_name += "Cn" + Cn_num;
		}
		if (Uut_num > 0) {
			final_name += "Uut" + Uut_num;
		}
		if (Fl_num > 0) {
			final_name += "Fl" + Fl_num;
		}
		if (Uup_num > 0) {
			final_name += "Uup" + Uup_num;
		}
		if (Lv_num > 0) {
			final_name += "Lv" + Lv_num;
		}
		if (Uus_num > 0) {
			final_name += "Uus" + Uus_num;
		}
		if (Uuo_num > 0) {
			final_name += "Uuo" + Uuo_num;
		}
		if (Xe_num > 0) {
			final_name += "Xe" + Xe_num;
		}
		return final_name;		
	}
	/**
	 * This function only works for C, H, N, O, P, S, F, Cl, Br formulas
	 * @param formula
	 * @return
	 */
	public static String standardize_name_DECOY(String formula) {
		if (!isNumeric(formula.substring(formula.length() - 1, formula.length()))) {
			formula = formula + "1";			
		}
		int Ci_num = retrieve_num_element(formula, "Ci");
		int H_num = retrieve_num_element(formula, "H") + 1;
		int He_num = retrieve_num_element(formula, "He");
		int Li_num = retrieve_num_element(formula, "Li");
		int Be_num = retrieve_num_element(formula, "Be");
		int B_num = retrieve_num_element(formula, "B");
		int C_num = retrieve_num_element(formula, "C");
		int N_num = retrieve_num_element(formula, "N");
		int O_num = retrieve_num_element(formula, "O");
		int F_num = retrieve_num_element(formula, "F");
		int Ne_num = retrieve_num_element(formula, "Ne");
		int Na_num = retrieve_num_element(formula, "Na");
		int Mg_num = retrieve_num_element(formula, "Mg");
		int Al_num = retrieve_num_element(formula, "Al");
		int Si_num = retrieve_num_element(formula, "Si");
		int P_num = retrieve_num_element(formula, "P");
		int S_num = retrieve_num_element(formula, "S");
		int Cl_num = retrieve_num_element(formula, "Cl");
		int Ar_num = retrieve_num_element(formula, "Ar");
		int K_num = retrieve_num_element(formula, "K");
		int Ca_num = retrieve_num_element(formula, "Ca");
		int Sc_num = retrieve_num_element(formula, "Sc");
		int Ti_num = retrieve_num_element(formula, "Ti");
		int V_num = retrieve_num_element(formula, "V");
		int Cr_num = retrieve_num_element(formula, "Cr");
		int Mn_num = retrieve_num_element(formula, "Mn");
		int Fe_num = retrieve_num_element(formula, "Fe");
		int Co_num = retrieve_num_element(formula, "Co");
		int Ni_num = retrieve_num_element(formula, "Ni");
		int Cu_num = retrieve_num_element(formula, "Cu");
		int Zn_num = retrieve_num_element(formula, "Zn");
		int Ga_num = retrieve_num_element(formula, "Ga");
		int Ge_num = retrieve_num_element(formula, "Ge");
		int As_num = retrieve_num_element(formula, "As");
		int Se_num = retrieve_num_element(formula, "Se");
		int Br_num = retrieve_num_element(formula, "Br");
		int Kr_num = retrieve_num_element(formula, "Kr");
		int Rb_num = retrieve_num_element(formula, "Rb");
		int Sr_num = retrieve_num_element(formula, "Sr");
		int Y_num = retrieve_num_element(formula, "Y");
		int Zr_num = retrieve_num_element(formula, "Zr");
		int Nb_num = retrieve_num_element(formula, "Nb");
		int Mo_num = retrieve_num_element(formula, "Mo");
		int Tc_num = retrieve_num_element(formula, "Tc");
		int Ru_num = retrieve_num_element(formula, "Ru");
		int Rh_num = retrieve_num_element(formula, "Rh");
		int Pd_num = retrieve_num_element(formula, "Pd");
		int Ag_num = retrieve_num_element(formula, "Ag");
		int Cd_num = retrieve_num_element(formula, "Cd");
		int In_num = retrieve_num_element(formula, "In");
		int Sn_num = retrieve_num_element(formula, "Sn");
		int Sb_num = retrieve_num_element(formula, "Sb");
		int Te_num = retrieve_num_element(formula, "Te");
		int I_num = retrieve_num_element(formula, "I");
		int Xe_num = retrieve_num_element(formula, "Xe");
		int Cs_num = retrieve_num_element(formula, "Cs");
		int Ba_num = retrieve_num_element(formula, "Ba");
		int La_num = retrieve_num_element(formula, "La");
		int Ce_num = retrieve_num_element(formula, "Ce");
		int Pr_num = retrieve_num_element(formula, "Pr");
		int Nd_num = retrieve_num_element(formula, "Nd");
		int Pm_num = retrieve_num_element(formula, "Pm");
		int Sm_num = retrieve_num_element(formula, "Sm");
		int Eu_num = retrieve_num_element(formula, "Eu");
		int Gd_num = retrieve_num_element(formula, "Gd");
		int Tb_num = retrieve_num_element(formula, "Tb");
		int Dy_num = retrieve_num_element(formula, "Dy");
		int Ho_num = retrieve_num_element(formula, "Ho");
		int Er_num = retrieve_num_element(formula, "Er");
		int Tm_num = retrieve_num_element(formula, "Tm");
		int Yb_num = retrieve_num_element(formula, "Yb");
		int Lu_num = retrieve_num_element(formula, "Lu");
		int Hf_num = retrieve_num_element(formula, "Hf");
		int Ta_num = retrieve_num_element(formula, "Ta");
		int W_num = retrieve_num_element(formula, "W");
		int Re_num = retrieve_num_element(formula, "Re");
		int Os_num = retrieve_num_element(formula, "Os");
		int Ir_num = retrieve_num_element(formula, "Ir");
		int Pt_num = retrieve_num_element(formula, "Pt");
		int Au_num = retrieve_num_element(formula, "Au");
		int Hg_num = retrieve_num_element(formula, "Hg");
		int Tl_num = retrieve_num_element(formula, "Tl");
		int Pb_num = retrieve_num_element(formula, "Pb");
		int Bi_num = retrieve_num_element(formula, "Bi");
		int Po_num = retrieve_num_element(formula, "Po");
		int At_num = retrieve_num_element(formula, "At");
		int Rn_num = retrieve_num_element(formula, "Rn");
		int Fr_num = retrieve_num_element(formula, "Fr");
		int Ra_num = retrieve_num_element(formula, "Ra");
		int Ac_num = retrieve_num_element(formula, "Ac");
		int Th_num = retrieve_num_element(formula, "Th");
		int Pa_num = retrieve_num_element(formula, "Pa");
		int U_num = retrieve_num_element(formula, "U");
		int Np_num = retrieve_num_element(formula, "Np");
		int Pu_num = retrieve_num_element(formula, "Pu");
		int Am_num = retrieve_num_element(formula, "Am");
		int Cm_num = retrieve_num_element(formula, "Cm");
		int Bk_num = retrieve_num_element(formula, "Bk");
		int Cf_num = retrieve_num_element(formula, "Cf");
		int Es_num = retrieve_num_element(formula, "Es");
		int Fm_num = retrieve_num_element(formula, "Fm");
		int Md_num = retrieve_num_element(formula, "Md");
		int No_num = retrieve_num_element(formula, "No");
		int Lr_num = retrieve_num_element(formula, "Lr");
		int Rf_num = retrieve_num_element(formula, "Rf");
		int Db_num = retrieve_num_element(formula, "Db");
		int Sg_num = retrieve_num_element(formula, "Sg");
		int Bh_num = retrieve_num_element(formula, "Bh");
		int Hs_num = retrieve_num_element(formula, "Hs");
		int Mt_num = retrieve_num_element(formula, "Mt");
		int Ds_num = retrieve_num_element(formula, "Ds");
		int Rg_num = retrieve_num_element(formula, "Rg");
		int Cn_num = retrieve_num_element(formula, "Cn");
		int Uut_num = retrieve_num_element(formula, "Uut");
		int Fl_num = retrieve_num_element(formula, "Fl");
		int Uup_num = retrieve_num_element(formula, "Uup");
		int Lv_num = retrieve_num_element(formula, "Lv");
		int Uus_num = retrieve_num_element(formula, "Uus");
		int Uuo_num = retrieve_num_element(formula, "Uuo");
		String final_name = "";
		if (Ci_num > 0) {
			final_name += "Ci" + Ci_num;
		}
		if (H_num > 0) {
			final_name += "H" + H_num;
		}
		if (He_num > 0) {
			final_name += "He" + He_num;
		}
		if (Li_num > 0) {
			final_name += "Li" + Li_num;
		}
		if (Be_num > 0) {
			final_name += "Be" + Be_num;
		}
		if (B_num > 0) {
			final_name += "B" + B_num;
		}
		if (C_num > 0) {
			final_name += "C" + C_num;
		}
		if (N_num > 0) {
			final_name += "N" + N_num;
		}
		if (O_num > 0) {
			final_name += "O" + O_num;
		}
		if (F_num > 0) {
			final_name += "F" + F_num;
		}
		if (Ne_num > 0) {
			final_name += "Ne" + Ne_num;
		}
		if (Na_num > 0) {
			final_name += "Na" + Na_num;
		}
		if (Mg_num > 0) {
			final_name += "Mg" + Mg_num;
		}
		if (Al_num > 0) {
			final_name += "Al" + Al_num;
		}
		if (Si_num > 0) {
			final_name += "Si" + Si_num;
		}
		if (P_num > 0) {
			final_name += "P" + P_num;
		}
		if (S_num > 0) {
			final_name += "S" + S_num;
		}
		if (Cl_num > 0) {
			final_name += "Cl" + Cl_num;
		}
		if (Ar_num > 0) {
			final_name += "Ar" + Ar_num;
		}
		if (K_num > 0) {
			final_name += "K" + K_num;
		}
		if (Ca_num > 0) {
			final_name += "Ca" + Ca_num;
		}
		if (Sc_num > 0) {
			final_name += "Sc" + Sc_num;
		}
		if (Ti_num > 0) {
			final_name += "Ti" + Ti_num;
		}
		if (V_num > 0) {
			final_name += "V" + V_num;
		}
		if (Cr_num > 0) {
			final_name += "Cr" + Cr_num;
		}
		if (Mn_num > 0) {
			final_name += "Mn" + Mn_num;
		}
		if (Fe_num > 0) {
			final_name += "Fe" + Fe_num;
		}
		if (Co_num > 0) {
			final_name += "Co" + Co_num;
		}
		if (Ni_num > 0) {
			final_name += "Ni" + Ni_num;
		}
		if (Cu_num > 0) {
			final_name += "Cu" + Cu_num;
		}
		if (Zn_num > 0) {
			final_name += "Zn" + Zn_num;
		}
		if (Ga_num > 0) {
			final_name += "Ga" + Ga_num;
		}
		if (Ge_num > 0) {
			final_name += "Ge" + Ge_num;
		}
		if (As_num > 0) {
			final_name += "As" + As_num;
		}
		if (Se_num > 0) {
			final_name += "Se" + Se_num;
		}
		if (Br_num > 0) {
			final_name += "Br" + Br_num;
		}
		if (Kr_num > 0) {
			final_name += "Kr" + Kr_num;
		}
		if (Rb_num > 0) {
			final_name += "Rb" + Rb_num;
		}
		if (Sr_num > 0) {
			final_name += "Sr" + Sr_num;
		}
		if (Y_num > 0) {
			final_name += "Y" + Y_num;
		}
		if (Zr_num > 0) {
			final_name += "Zr" + Zr_num;
		}
		if (Nb_num > 0) {
			final_name += "Nb" + Nb_num;
		}
		if (Mo_num > 0) {
			final_name += "Mo" + Mo_num;
		}
		if (Tc_num > 0) {
			final_name += "Tc" + Tc_num;
		}
		if (Ru_num > 0) {
			final_name += "Ru" + Ru_num;
		}
		if (Rh_num > 0) {
			final_name += "Rh" + Rh_num;
		}
		if (Pd_num > 0) {
			final_name += "Pd" + Pd_num;
		}
		if (Ag_num > 0) {
			final_name += "Ag" + Ag_num;
		}
		if (Cd_num > 0) {
			final_name += "Cd" + Cd_num;
		}
		if (In_num > 0) {
			final_name += "In" + In_num;
		}
		if (Sn_num > 0) {
			final_name += "Sn" + Sn_num;
		}
		if (Sb_num > 0) {
			final_name += "Sb" + Sb_num;
		}
		if (Te_num > 0) {
			final_name += "Te" + Te_num;
		}
		if (I_num > 0) {
			final_name += "I" + I_num;
		}
		if (Xe_num > 0) {
			final_name += "Xe" + Xe_num;
		}
		if (Cs_num > 0) {
			final_name += "Cs" + Cs_num;
		}
		if (Ba_num > 0) {
			final_name += "Ba" + Ba_num;
		}
		if (La_num > 0) {
			final_name += "La" + La_num;
		}
		if (Ce_num > 0) {
			final_name += "Ce" + Ce_num;
		}
		if (Pr_num > 0) {
			final_name += "Pr" + Pr_num;
		}
		if (Nd_num > 0) {
			final_name += "Nd" + Nd_num;
		}
		if (Pm_num > 0) {
			final_name += "Pm" + Pm_num;
		}
		if (Sm_num > 0) {
			final_name += "Sm" + Sm_num;
		}
		if (Eu_num > 0) {
			final_name += "Eu" + Eu_num;
		}
		if (Gd_num > 0) {
			final_name += "Gd" + Gd_num;
		}
		if (Tb_num > 0) {
			final_name += "Tb" + Tb_num;
		}
		if (Dy_num > 0) {
			final_name += "Dy" + Dy_num;
		}
		if (Ho_num > 0) {
			final_name += "Ho" + Ho_num;
		}
		if (Er_num > 0) {
			final_name += "Er" + Er_num;
		}
		if (Tm_num > 0) {
			final_name += "Tm" + Tm_num;
		}
		if (Yb_num > 0) {
			final_name += "Yb" + Yb_num;
		}
		if (Lu_num > 0) {
			final_name += "Lu" + Lu_num;
		}
		if (Hf_num > 0) {
			final_name += "Hf" + Hf_num;
		}
		if (Ta_num > 0) {
			final_name += "Ta" + Ta_num;
		}
		if (W_num > 0) {
			final_name += "W" + W_num;
		}
		if (Re_num > 0) {
			final_name += "Re" + Re_num;
		}
		if (Os_num > 0) {
			final_name += "Os" + Os_num;
		}
		if (Ir_num > 0) {
			final_name += "Ir" + Ir_num;
		}
		if (Pt_num > 0) {
			final_name += "Pt" + Pt_num;
		}
		if (Au_num > 0) {
			final_name += "Au" + Au_num;
		}
		if (Hg_num > 0) {
			final_name += "Hg" + Hg_num;
		}
		if (Tl_num > 0) {
			final_name += "Tl" + Tl_num;
		}
		if (Pb_num > 0) {
			final_name += "Pb" + Pb_num;
		}
		if (Bi_num > 0) {
			final_name += "Bi" + Bi_num;
		}
		if (Po_num > 0) {
			final_name += "Po" + Po_num;
		}
		if (At_num > 0) {
			final_name += "At" + At_num;
		}
		if (Rn_num > 0) {
			final_name += "Rn" + Rn_num;
		}
		if (Fr_num > 0) {
			final_name += "Fr" + Fr_num;
		}
		if (Ra_num > 0) {
			final_name += "Ra" + Ra_num;
		}
		if (Ac_num > 0) {
			final_name += "Ac" + Ac_num;
		}
		if (Th_num > 0) {
			final_name += "Th" + Th_num;
		}
		if (Pa_num > 0) {
			final_name += "Pa" + Pa_num;
		}
		if (U_num > 0) {
			final_name += "U" + U_num;
		}
		if (Np_num > 0) {
			final_name += "Np" + Np_num;
		}
		if (Pu_num > 0) {
			final_name += "Pu" + Pu_num;
		}
		if (Am_num > 0) {
			final_name += "Am" + Am_num;
		}
		if (Cm_num > 0) {
			final_name += "Cm" + Cm_num;
		}
		if (Bk_num > 0) {
			final_name += "Bk" + Bk_num;
		}
		if (Cf_num > 0) {
			final_name += "Cf" + Cf_num;
		}
		if (Es_num > 0) {
			final_name += "Es" + Es_num;
		}
		if (Fm_num > 0) {
			final_name += "Fm" + Fm_num;
		}
		if (Md_num > 0) {
			final_name += "Md" + Md_num;
		}
		if (No_num > 0) {
			final_name += "No" + No_num;
		}
		if (Lr_num > 0) {
			final_name += "Lr" + Lr_num;
		}
		if (Rf_num > 0) {
			final_name += "Rf" + Rf_num;
		}
		if (Db_num > 0) {
			final_name += "Db" + Db_num;
		}
		if (Sg_num > 0) {
			final_name += "Sg" + Sg_num;
		}
		if (Bh_num > 0) {
			final_name += "Bh" + Bh_num;
		}
		if (Hs_num > 0) {
			final_name += "Hs" + Hs_num;
		}
		if (Mt_num > 0) {
			final_name += "Mt" + Mt_num;
		}
		if (Ds_num > 0) {
			final_name += "Ds" + Ds_num;
		}
		if (Rg_num > 0) {
			final_name += "Rg" + Rg_num;
		}
		if (Cn_num > 0) {
			final_name += "Cn" + Cn_num;
		}
		if (Uut_num > 0) {
			final_name += "Uut" + Uut_num;
		}
		if (Fl_num > 0) {
			final_name += "Fl" + Fl_num;
		}
		if (Uup_num > 0) {
			final_name += "Uup" + Uup_num;
		}
		if (Lv_num > 0) {
			final_name += "Lv" + Lv_num;
		}
		if (Uus_num > 0) {
			final_name += "Uus" + Uus_num;
		}
		if (Uuo_num > 0) {
			final_name += "Uuo" + Uuo_num;
		}
		return final_name;
	}
	/**
	 * This function only works for C, H, N, O, P, S, F, Cl, Br formulas
	 * @param formula
	 * @return
	 */
	public static boolean check_formula_valid_element_corrected(String formula) {
		
		formula = standardize_name(formula);		
		if (formula.trim().equals("")) {
			return false;
		}
		int Ci_num = retrieve_num_element(formula, "Ci");
		int H_num = retrieve_num_element(formula, "H");
		int He_num = retrieve_num_element(formula, "He");
		int Li_num = retrieve_num_element(formula, "Li");
		int Be_num = retrieve_num_element(formula, "Be");
		int B_num = retrieve_num_element(formula, "B");
		int C_num = retrieve_num_element(formula, "C");
		int N_num = retrieve_num_element(formula, "N");
		int O_num = retrieve_num_element(formula, "O");
		int F_num = retrieve_num_element(formula, "F");
		int Ne_num = retrieve_num_element(formula, "Ne");
		int Na_num = retrieve_num_element(formula, "Na");
		int Mg_num = retrieve_num_element(formula, "Mg");
		int Al_num = retrieve_num_element(formula, "Al");
		int Si_num = retrieve_num_element(formula, "Si");
		int P_num = retrieve_num_element(formula, "P");
		int S_num = retrieve_num_element(formula, "S");
		int Cl_num = retrieve_num_element(formula, "Cl");
		int Ar_num = retrieve_num_element(formula, "Ar");
		int K_num = retrieve_num_element(formula, "K");
		int Ca_num = retrieve_num_element(formula, "Ca");
		int Sc_num = retrieve_num_element(formula, "Sc");
		int Ti_num = retrieve_num_element(formula, "Ti");
		int V_num = retrieve_num_element(formula, "V");
		int Cr_num = retrieve_num_element(formula, "Cr");
		int Mn_num = retrieve_num_element(formula, "Mn");
		int Fe_num = retrieve_num_element(formula, "Fe");
		int Co_num = retrieve_num_element(formula, "Co");
		int Ni_num = retrieve_num_element(formula, "Ni");
		int Cu_num = retrieve_num_element(formula, "Cu");
		int Zn_num = retrieve_num_element(formula, "Zn");
		int Ga_num = retrieve_num_element(formula, "Ga");
		int Ge_num = retrieve_num_element(formula, "Ge");
		int As_num = retrieve_num_element(formula, "As");
		int Se_num = retrieve_num_element(formula, "Se");
		int Br_num = retrieve_num_element(formula, "Br");
		int Kr_num = retrieve_num_element(formula, "Kr");
		int Rb_num = retrieve_num_element(formula, "Rb");
		int Sr_num = retrieve_num_element(formula, "Sr");
		int Y_num = retrieve_num_element(formula, "Y");
		int Zr_num = retrieve_num_element(formula, "Zr");
		int Nb_num = retrieve_num_element(formula, "Nb");
		int Mo_num = retrieve_num_element(formula, "Mo");
		int Tc_num = retrieve_num_element(formula, "Tc");
		int Ru_num = retrieve_num_element(formula, "Ru");
		int Rh_num = retrieve_num_element(formula, "Rh");
		int Pd_num = retrieve_num_element(formula, "Pd");
		int Ag_num = retrieve_num_element(formula, "Ag");
		int Cd_num = retrieve_num_element(formula, "Cd");
		int In_num = retrieve_num_element(formula, "In");
		int Sn_num = retrieve_num_element(formula, "Sn");
		int Sb_num = retrieve_num_element(formula, "Sb");
		int Te_num = retrieve_num_element(formula, "Te");
		int I_num = retrieve_num_element(formula, "I");
		int Xe_num = retrieve_num_element(formula, "Xe");
		int Cs_num = retrieve_num_element(formula, "Cs");
		int Ba_num = retrieve_num_element(formula, "Ba");
		int La_num = retrieve_num_element(formula, "La");
		int Ce_num = retrieve_num_element(formula, "Ce");
		int Pr_num = retrieve_num_element(formula, "Pr");
		int Nd_num = retrieve_num_element(formula, "Nd");
		int Pm_num = retrieve_num_element(formula, "Pm");
		int Sm_num = retrieve_num_element(formula, "Sm");
		int Eu_num = retrieve_num_element(formula, "Eu");
		int Gd_num = retrieve_num_element(formula, "Gd");
		int Tb_num = retrieve_num_element(formula, "Tb");
		int Dy_num = retrieve_num_element(formula, "Dy");
		int Ho_num = retrieve_num_element(formula, "Ho");
		int Er_num = retrieve_num_element(formula, "Er");
		int Tm_num = retrieve_num_element(formula, "Tm");
		int Yb_num = retrieve_num_element(formula, "Yb");
		int Lu_num = retrieve_num_element(formula, "Lu");
		int Hf_num = retrieve_num_element(formula, "Hf");
		int Ta_num = retrieve_num_element(formula, "Ta");
		int W_num = retrieve_num_element(formula, "W");
		int Re_num = retrieve_num_element(formula, "Re");
		int Os_num = retrieve_num_element(formula, "Os");
		int Ir_num = retrieve_num_element(formula, "Ir");
		int Pt_num = retrieve_num_element(formula, "Pt");
		int Au_num = retrieve_num_element(formula, "Au");
		int Hg_num = retrieve_num_element(formula, "Hg");
		int Tl_num = retrieve_num_element(formula, "Tl");
		int Pb_num = retrieve_num_element(formula, "Pb");
		int Bi_num = retrieve_num_element(formula, "Bi");
		int Po_num = retrieve_num_element(formula, "Po");
		int At_num = retrieve_num_element(formula, "At");
		int Rn_num = retrieve_num_element(formula, "Rn");
		int Fr_num = retrieve_num_element(formula, "Fr");
		int Ra_num = retrieve_num_element(formula, "Ra");
		int Ac_num = retrieve_num_element(formula, "Ac");
		int Th_num = retrieve_num_element(formula, "Th");
		int Pa_num = retrieve_num_element(formula, "Pa");
		int U_num = retrieve_num_element(formula, "U");
		int Np_num = retrieve_num_element(formula, "Np");
		int Pu_num = retrieve_num_element(formula, "Pu");
		int Am_num = retrieve_num_element(formula, "Am");
		int Cm_num = retrieve_num_element(formula, "Cm");
		int Bk_num = retrieve_num_element(formula, "Bk");
		int Cf_num = retrieve_num_element(formula, "Cf");
		int Es_num = retrieve_num_element(formula, "Es");
		int Fm_num = retrieve_num_element(formula, "Fm");
		int Md_num = retrieve_num_element(formula, "Md");
		int No_num = retrieve_num_element(formula, "No");
		int Lr_num = retrieve_num_element(formula, "Lr");
		int Rf_num = retrieve_num_element(formula, "Rf");
		int Db_num = retrieve_num_element(formula, "Db");
		int Sg_num = retrieve_num_element(formula, "Sg");
		int Bh_num = retrieve_num_element(formula, "Bh");
		int Hs_num = retrieve_num_element(formula, "Hs");
		int Mt_num = retrieve_num_element(formula, "Mt");
		int Ds_num = retrieve_num_element(formula, "Ds");
		int Rg_num = retrieve_num_element(formula, "Rg");
		int Cn_num = retrieve_num_element(formula, "Cn");
		int Uut_num = retrieve_num_element(formula, "Uut");
		int Fl_num = retrieve_num_element(formula, "Fl");
		int Uup_num = retrieve_num_element(formula, "Uup");
		int Lv_num = retrieve_num_element(formula, "Lv");
		int Uus_num = retrieve_num_element(formula, "Uus");
		int Uuo_num = retrieve_num_element(formula, "Uuo");
		String final_name = "";
		
		int total = 0;
		if (Ci_num > 0) {
			final_name += "Ci" + Ci_num;
			if (Ci_num > 0) {
				return false;
			}
		}
		if (H_num > 0) {
			final_name += "H" + H_num;
			
		}
		if (He_num > 0) {
			final_name += "He" + He_num;
			if (He_num > 0) {
				return false;
			}
		}
		if (Li_num > 0) {
			final_name += "Li" + Li_num;
			if (Li_num > 0) {
				return false;
			}
		}
		if (Be_num > 0) {
			final_name += "Be" + Be_num;
			if (Be_num > 0) {
				return false;
			}
		}
		if (B_num > 0) {
			final_name += "B" + B_num;
			
		}
		if (C_num > 0) {
			final_name += "C" + C_num;
			
		}
		if (N_num > 0) {
			final_name += "N" + N_num;
			
		}
		if (O_num > 0) {
			final_name += "O" + O_num;
			
		}
		if (F_num > 0) {
			final_name += "F" + F_num;
			
		}
		if (Ne_num > 0) {
			final_name += "Ne" + Ne_num;
			if (Ne_num > 0) {
				return false;
			}
		}
		if (Na_num > 0) {
			final_name += "Na" + Na_num;
			
		}
		if (Mg_num > 0) {
			final_name += "Mg" + Mg_num;
			
		}
		if (Al_num > 0) {
			final_name += "Al" + Al_num;
			
		}
		if (Si_num > 0) {
			final_name += "Si" + Si_num;
			
		}
		if (P_num > 0) {
			final_name += "P" + P_num;
			if (P_num > 300) {
				return false;
			}
		}
		if (S_num > 0) {
			final_name += "S" + S_num;
			if (S_num > 300) {
				return false;
			}
		}
		if (Cl_num > 0) {
			final_name += "Cl" + Cl_num;
		}
		if (Ar_num > 0) {
			final_name += "Ar" + Ar_num;
			if (Ar_num > 0) {
				return false;
			}
		}
		if (K_num > 0) {
			final_name += "K" + K_num;
			
		}
		if (Ca_num > 0) {
			final_name += "Ca" + Ca_num;
			
		}
		if (Sc_num > 0) {
			final_name += "Sc" + Sc_num;
			if (Sc_num > 0) {
				return false;
			}
		}
		if (Ti_num > 0) {
			final_name += "Ti" + Ti_num;
			if (Ti_num > 0) {
				return false;
			}
		}
		if (V_num > 0) {
			final_name += "V" + V_num;
			if (V_num > 0) {
				return false;
			}
		}
		if (Cr_num > 0) {
			final_name += "Cr" + Cr_num;
			if (Cr_num > 0) {
				return false;
			}
		}
		if (Mn_num > 0) {
			final_name += "Mn" + Mn_num;
			if (Mn_num > 0) {
				return false;
			}
		}
		if (Fe_num > 0) {
			final_name += "Fe" + Fe_num;
			if (Fe_num > 0) {
				return false;
			}
		}
		if (Co_num > 0) {
			final_name += "Co" + Co_num;
			if (Co_num > 0) {
				return false;
			}
		}
		if (Ni_num > 0) {
			final_name += "Ni" + Ni_num;
			if (Ni_num > 0) {
				return false;
			}
		}
		if (Cu_num > 0) {
			final_name += "Cu" + Cu_num;
			if (Cu_num > 0) {
				return false;
			}
		}
		if (Zn_num > 0) {
			final_name += "Zn" + Zn_num;
			if (Zn_num > 0) {
				return false;
			}
		}
		if (Ga_num > 0) {
			final_name += "Ga" + Ga_num;
			if (Ga_num > 0) {
				return false;
			}
		}
		if (Ge_num > 0) {
			final_name += "Ge" + Ge_num;
			if (Ge_num > 0) {
				return false;
			}
		}
		if (As_num > 0) {
			final_name += "As" + As_num;
			if (As_num > 0) {
				return false;
			}
		}
		if (Se_num > 0) {
			final_name += "Se" + Se_num;
			if (Se_num > 0) {
				return false;
			}
		}
		if (Br_num > 0) {
			final_name += "Br" + Br_num;
		}
		if (Kr_num > 0) {
			final_name += "Kr" + Kr_num;
			if (Kr_num > 0) {
				return false;
			}
		}
		if (Rb_num > 0) {
			final_name += "Rb" + Rb_num;
			if (Rb_num > 0) {
				return false;
			}
		}
		if (Sr_num > 0) {
			final_name += "Sr" + Sr_num;
			if (Sr_num > 0) {
				return false;
			}
		}
		if (Y_num > 0) {
			final_name += "Y" + Y_num;
			if (Y_num > 0) {
				return false;
			}
		}
		if (Zr_num > 0) {
			final_name += "Zr" + Zr_num;
			if (Zr_num > 0) {
				return false;
			}
		}
		if (Nb_num > 0) {
			final_name += "Nb" + Nb_num;
			if (Nb_num > 0) {
				return false;
			}
		}
		if (Mo_num > 0) {
			final_name += "Mo" + Mo_num;
			if (Mo_num > 0) {
				return false;
			}
		}
		if (Tc_num > 0) {
			final_name += "Tc" + Tc_num;
			if (Tc_num > 0) {
				return false;
			}
		}
		if (Ru_num > 0) {
			final_name += "Ru" + Ru_num;
			if (Ru_num > 0) {
				return false;
			}
		}
		if (Rh_num > 0) {
			final_name += "Rh" + Rh_num;
			if (Rh_num > 0) {
				return false;
			}
		}
		if (Pd_num > 0) {
			final_name += "Pd" + Pd_num;
			if (Pd_num > 0) {
				return false;
			}
		}
		if (Ag_num > 0) {
			final_name += "Ag" + Ag_num;
			if (Ag_num > 0) {
				return false;
			}
		}
		if (Cd_num > 0) {
			final_name += "Cd" + Cd_num;
			if (Cd_num > 0) {
				return false;
			}
		}
		if (In_num > 0) {
			final_name += "In" + In_num;
			if (In_num > 0) {
				return false;
			}
		}
		if (Sn_num > 0) {
			final_name += "Sn" + Sn_num;
			if (Sn_num > 0) {
				return false;
			}
		}
		if (Sb_num > 0) {
			final_name += "Sb" + Sb_num;
			if (Sb_num > 0) {
				return false;
			}
		}
		if (Te_num > 0) {
			final_name += "Te" + Te_num;
			if (Te_num > 0) {
				return false;
			}
		}
		if (I_num > 0) {
			final_name += "I" + I_num;
			
		}
		if (Xe_num > 0) {
			final_name += "Xe" + Xe_num;
			if (Xe_num > 0) {
				return false;
			}
		}
		if (Cs_num > 0) {
			final_name += "Cs" + Cs_num;
			if (Cs_num > 0) {
				return false;
			}
		}
		if (Ba_num > 0) {
			final_name += "Ba" + Ba_num;
			if (Ba_num > 0) {
				return false;
			}
		}
		if (La_num > 0) {
			final_name += "La" + La_num;
			if (La_num > 0) {
				return false;
			}
		}
		if (Ce_num > 0) {
			final_name += "Ce" + Ce_num;
			if (Ce_num > 0) {
				return false;
			}
		}
		if (Pr_num > 0) {
			final_name += "Pr" + Pr_num;
			if (Pr_num > 0) {
				return false;
			}
		}
		if (Nd_num > 0) {
			final_name += "Nd" + Nd_num;
			if (Nd_num > 0) {
				return false;
			}
		}
		if (Pm_num > 0) {
			final_name += "Pm" + Pm_num;
			if (Pm_num > 0) {
				return false;
			}
		}
		if (Sm_num > 0) {
			final_name += "Sm" + Sm_num;
			if (Sm_num > 0) {
				return false;
			}
		}
		if (Eu_num > 0) {
			final_name += "Eu" + Eu_num;
			if (Eu_num > 0) {
				return false;
			}
		}
		if (Gd_num > 0) {
			final_name += "Gd" + Gd_num;
			if (Gd_num > 0) {
				return false;
			}
		}
		if (Tb_num > 0) {
			final_name += "Tb" + Tb_num;
			if (Tb_num > 0) {
				return false;
			}
		}
		if (Dy_num > 0) {
			final_name += "Dy" + Dy_num;
			if (Dy_num > 0) {
				return false;
			}
		}
		if (Ho_num > 0) {
			final_name += "Ho" + Ho_num;
			if (Ho_num > 0) {
				return false;
			}
		}
		if (Er_num > 0) {
			final_name += "Er" + Er_num;
			if (Er_num > 0) {
				return false;
			}
		}
		if (Tm_num > 0) {
			final_name += "Tm" + Tm_num;
			if (Tm_num > 0) {
				return false;
			}
		}
		if (Yb_num > 0) {
			final_name += "Yb" + Yb_num;
			if (Yb_num > 0) {
				return false;
			}
		}
		if (Lu_num > 0) {
			final_name += "Lu" + Lu_num;
			if (Lu_num > 0) {
				return false;
			}
		}
		if (Hf_num > 0) {
			final_name += "Hf" + Hf_num;
			if (Hf_num > 0) {
				return false;
			}
		}
		if (Ta_num > 0) {
			final_name += "Ta" + Ta_num;
			if (Ta_num > 0) {
				return false;
			}
		}
		if (W_num > 0) {
			final_name += "W" + W_num;
			if (W_num > 0) {
				return false;
			}
		}
		if (Re_num > 0) {
			final_name += "Re" + Re_num;
			if (Re_num > 0) {
				return false;
			}
		}
		if (Os_num > 0) {
			final_name += "Os" + Os_num;
			if (Os_num > 0) {
				return false;
			}
		}
		if (Ir_num > 0) {
			final_name += "Ir" + Ir_num;
			if (Ir_num > 0) {
				return false;
			}
		}
		if (Pt_num > 0) {
			final_name += "Pt" + Pt_num;
			if (Pt_num > 0) {
				return false;
			}
		}
		if (Au_num > 0) {
			final_name += "Au" + Au_num;
			if (Au_num > 0) {
				return false;
			}
		}
		if (Hg_num > 0) {
			final_name += "Hg" + Hg_num;
			if (Hg_num > 0) {
				return false;
			}
		}
		if (Tl_num > 0) {
			final_name += "Tl" + Tl_num;
			if (Tl_num > 0) {
				return false;
			}
		}
		if (Pb_num > 0) {
			final_name += "Pb" + Pb_num;
			if (Pb_num > 0) {
				return false;
			}
		}
		if (Bi_num > 0) {
			final_name += "Bi" + Bi_num;
			if (Bi_num > 0) {
				return false;
			}
		}
		if (Po_num > 0) {
			final_name += "Po" + Po_num;
			if (Po_num > 0) {
				return false;
			}
		}
		if (At_num > 0) {
			final_name += "At" + At_num;
			if (At_num > 0) {
				return false;
			}
		}
		if (Rn_num > 0) {
			final_name += "Rn" + Rn_num;
			if (Rn_num > 0) {
				return false;
			}
		}
		if (Fr_num > 0) {
			final_name += "Fr" + Fr_num;
			if (Fr_num > 0) {
				return false;
			}
		}
		if (Ra_num > 0) {
			final_name += "Ra" + Ra_num;
			if (Ra_num > 0) {
				return false;
			}
		}
		if (Ac_num > 0) {
			final_name += "Ac" + Ac_num;
			if (Ac_num > 0) {
				return false;
			}
		}
		if (Th_num > 0) {
			final_name += "Th" + Th_num;
			if (Th_num > 0) {
				return false;
			}
		}
		if (Pa_num > 0) {
			final_name += "Pa" + Pa_num;
			if (Pa_num > 0) {
				return false;
			}
		}
		if (U_num > 0) {
			final_name += "U" + U_num;
			if (U_num > 0) {
				return false;
			}
		}
		if (Np_num > 0) {
			final_name += "Np" + Np_num;
			if (Np_num > 0) {
				return false;
			}
		}
		if (Pu_num > 0) {
			final_name += "Pu" + Pu_num;
			if (Pu_num > 0) {
				return false;
			}
		}
		if (Am_num > 0) {
			final_name += "Am" + Am_num;
			if (Am_num > 0) {
				return false;
			}
		}
		if (Cm_num > 0) {
			final_name += "Cm" + Cm_num;
			if (Cm_num > 0) {
				return false;
			}
		}
		if (Bk_num > 0) {
			final_name += "Bk" + Bk_num;
			if (Bk_num > 0) {
				return false;
			}
		}
		if (Cf_num > 0) {
			final_name += "Cf" + Cf_num;
			if (Cf_num > 0) {
				return false;
			}
		}
		if (Es_num > 0) {
			final_name += "Es" + Es_num;
			if (Es_num > 0) {
				return false;
			}
		}
		if (Fm_num > 0) {
			final_name += "Fm" + Fm_num;
			if (Fm_num > 0) {
				return false;
			}
		}
		if (Md_num > 0) {
			final_name += "Md" + Md_num;
			if (Md_num > 0) {
				return false;
			}
		}
		if (No_num > 0) {
			final_name += "No" + No_num;
			if (No_num > 0) {
				return false;
			}
		}
		if (Lr_num > 0) {
			final_name += "Lr" + Lr_num;
			if (Lr_num > 0) {
				return false;
			}
		}
		if (Rf_num > 0) {
			final_name += "Rf" + Rf_num;
			if (Rf_num > 0) {
				return false;
			}
		}
		if (Db_num > 0) {
			final_name += "Db" + Db_num;
			if (Db_num > 0) {
				return false;
			}
		}
		if (Sg_num > 0) {
			final_name += "Sg" + Sg_num;
			if (Sg_num > 0) {
				return false;
			}
		}
		if (Bh_num > 0) {
			final_name += "Bh" + Bh_num;
			if (Bh_num > 0) {
				return false;
			}
		}
		if (Hs_num > 0) {
			final_name += "Hs" + Hs_num;
			if (Hs_num > 0) {
				return false;
			}
		}
		if (Mt_num > 0) {
			final_name += "Mt" + Mt_num;
			if (Mt_num > 0) {
				return false;
			}
		}
		if (Ds_num > 0) {
			final_name += "Ds" + Ds_num;
			if (Ds_num > 0) {
				return false;
			}
		}
		if (Rg_num > 0) {
			final_name += "Rg" + Rg_num;
			if (Rg_num > 0) {
				return false;
			}
		}
		if (Cn_num > 0) {
			final_name += "Cn" + Cn_num;
			if (Cn_num > 0) {
				return false;
			}
		}
		if (Uut_num > 0) {
			final_name += "Uut" + Uut_num;
			if (Uut_num > 0) {
				return false;
			}
		}
		if (Fl_num > 0) {
			final_name += "Fl" + Fl_num;
			if (Fl_num > 0) {
				return false;
			}
		}
		if (Uup_num > 0) {
			final_name += "Uup" + Uup_num;
			if (Uup_num > 0) {
				return false;
			}
		}
		if (Lv_num > 0) {
			final_name += "Lv" + Lv_num;
			if (Lv_num > 0) {
				return false;
			}
		}
		if (Uus_num > 0) {
			final_name += "Uus" + Uus_num;
			if (Uus_num > 0) {
				return false;
			}
		}
		if (Uuo_num > 0) {
			final_name += "Uuo" + Uuo_num;
			if (Uuo_num > 0) {
				return false;
			}
		}
		return true;
	}
	/**
	 * This function only works for C, H, N, O, P, S, F, Cl, Br formulas
	 * @param formula
	 * @return
	 */
	public static boolean check_formula_valid_element(String formula) {
		if (formula.trim().equals("")) {
			return false;
		}
		int Ci_num = retrieve_num_element(formula, "Ci");
		int H_num = retrieve_num_element(formula, "H");
		int He_num = retrieve_num_element(formula, "He");
		int Li_num = retrieve_num_element(formula, "Li");
		int Be_num = retrieve_num_element(formula, "Be");
		int B_num = retrieve_num_element(formula, "B");
		int C_num = retrieve_num_element(formula, "C");
		int N_num = retrieve_num_element(formula, "N");
		int O_num = retrieve_num_element(formula, "O");
		int F_num = retrieve_num_element(formula, "F");
		int Ne_num = retrieve_num_element(formula, "Ne");
		int Na_num = retrieve_num_element(formula, "Na");
		int Mg_num = retrieve_num_element(formula, "Mg");
		int Al_num = retrieve_num_element(formula, "Al");
		int Si_num = retrieve_num_element(formula, "Si");
		int P_num = retrieve_num_element(formula, "P");
		int S_num = retrieve_num_element(formula, "S");
		int Cl_num = retrieve_num_element(formula, "Cl");
		int Ar_num = retrieve_num_element(formula, "Ar");
		int K_num = retrieve_num_element(formula, "K");
		int Ca_num = retrieve_num_element(formula, "Ca");
		int Sc_num = retrieve_num_element(formula, "Sc");
		int Ti_num = retrieve_num_element(formula, "Ti");
		int V_num = retrieve_num_element(formula, "V");
		int Cr_num = retrieve_num_element(formula, "Cr");
		int Mn_num = retrieve_num_element(formula, "Mn");
		int Fe_num = retrieve_num_element(formula, "Fe");
		int Co_num = retrieve_num_element(formula, "Co");
		int Ni_num = retrieve_num_element(formula, "Ni");
		int Cu_num = retrieve_num_element(formula, "Cu");
		int Zn_num = retrieve_num_element(formula, "Zn");
		int Ga_num = retrieve_num_element(formula, "Ga");
		int Ge_num = retrieve_num_element(formula, "Ge");
		int As_num = retrieve_num_element(formula, "As");
		int Se_num = retrieve_num_element(formula, "Se");
		int Br_num = retrieve_num_element(formula, "Br");
		int Kr_num = retrieve_num_element(formula, "Kr");
		int Rb_num = retrieve_num_element(formula, "Rb");
		int Sr_num = retrieve_num_element(formula, "Sr");
		int Y_num = retrieve_num_element(formula, "Y");
		int Zr_num = retrieve_num_element(formula, "Zr");
		int Nb_num = retrieve_num_element(formula, "Nb");
		int Mo_num = retrieve_num_element(formula, "Mo");
		int Tc_num = retrieve_num_element(formula, "Tc");
		int Ru_num = retrieve_num_element(formula, "Ru");
		int Rh_num = retrieve_num_element(formula, "Rh");
		int Pd_num = retrieve_num_element(formula, "Pd");
		int Ag_num = retrieve_num_element(formula, "Ag");
		int Cd_num = retrieve_num_element(formula, "Cd");
		int In_num = retrieve_num_element(formula, "In");
		int Sn_num = retrieve_num_element(formula, "Sn");
		int Sb_num = retrieve_num_element(formula, "Sb");
		int Te_num = retrieve_num_element(formula, "Te");
		int I_num = retrieve_num_element(formula, "I");
		int Xe_num = retrieve_num_element(formula, "Xe");
		int Cs_num = retrieve_num_element(formula, "Cs");
		int Ba_num = retrieve_num_element(formula, "Ba");
		int La_num = retrieve_num_element(formula, "La");
		int Ce_num = retrieve_num_element(formula, "Ce");
		int Pr_num = retrieve_num_element(formula, "Pr");
		int Nd_num = retrieve_num_element(formula, "Nd");
		int Pm_num = retrieve_num_element(formula, "Pm");
		int Sm_num = retrieve_num_element(formula, "Sm");
		int Eu_num = retrieve_num_element(formula, "Eu");
		int Gd_num = retrieve_num_element(formula, "Gd");
		int Tb_num = retrieve_num_element(formula, "Tb");
		int Dy_num = retrieve_num_element(formula, "Dy");
		int Ho_num = retrieve_num_element(formula, "Ho");
		int Er_num = retrieve_num_element(formula, "Er");
		int Tm_num = retrieve_num_element(formula, "Tm");
		int Yb_num = retrieve_num_element(formula, "Yb");
		int Lu_num = retrieve_num_element(formula, "Lu");
		int Hf_num = retrieve_num_element(formula, "Hf");
		int Ta_num = retrieve_num_element(formula, "Ta");
		int W_num = retrieve_num_element(formula, "W");
		int Re_num = retrieve_num_element(formula, "Re");
		int Os_num = retrieve_num_element(formula, "Os");
		int Ir_num = retrieve_num_element(formula, "Ir");
		int Pt_num = retrieve_num_element(formula, "Pt");
		int Au_num = retrieve_num_element(formula, "Au");
		int Hg_num = retrieve_num_element(formula, "Hg");
		int Tl_num = retrieve_num_element(formula, "Tl");
		int Pb_num = retrieve_num_element(formula, "Pb");
		int Bi_num = retrieve_num_element(formula, "Bi");
		int Po_num = retrieve_num_element(formula, "Po");
		int At_num = retrieve_num_element(formula, "At");
		int Rn_num = retrieve_num_element(formula, "Rn");
		int Fr_num = retrieve_num_element(formula, "Fr");
		int Ra_num = retrieve_num_element(formula, "Ra");
		int Ac_num = retrieve_num_element(formula, "Ac");
		int Th_num = retrieve_num_element(formula, "Th");
		int Pa_num = retrieve_num_element(formula, "Pa");
		int U_num = retrieve_num_element(formula, "U");
		int Np_num = retrieve_num_element(formula, "Np");
		int Pu_num = retrieve_num_element(formula, "Pu");
		int Am_num = retrieve_num_element(formula, "Am");
		int Cm_num = retrieve_num_element(formula, "Cm");
		int Bk_num = retrieve_num_element(formula, "Bk");
		int Cf_num = retrieve_num_element(formula, "Cf");
		int Es_num = retrieve_num_element(formula, "Es");
		int Fm_num = retrieve_num_element(formula, "Fm");
		int Md_num = retrieve_num_element(formula, "Md");
		int No_num = retrieve_num_element(formula, "No");
		int Lr_num = retrieve_num_element(formula, "Lr");
		int Rf_num = retrieve_num_element(formula, "Rf");
		int Db_num = retrieve_num_element(formula, "Db");
		int Sg_num = retrieve_num_element(formula, "Sg");
		int Bh_num = retrieve_num_element(formula, "Bh");
		int Hs_num = retrieve_num_element(formula, "Hs");
		int Mt_num = retrieve_num_element(formula, "Mt");
		int Ds_num = retrieve_num_element(formula, "Ds");
		int Rg_num = retrieve_num_element(formula, "Rg");
		int Cn_num = retrieve_num_element(formula, "Cn");
		int Uut_num = retrieve_num_element(formula, "Uut");
		int Fl_num = retrieve_num_element(formula, "Fl");
		int Uup_num = retrieve_num_element(formula, "Uup");
		int Lv_num = retrieve_num_element(formula, "Lv");
		int Uus_num = retrieve_num_element(formula, "Uus");
		int Uuo_num = retrieve_num_element(formula, "Uuo");
		String final_name = "";
		
		int total = 0;
		if (Ci_num > 0) {
			final_name += "Ci" + Ci_num;			
		}
		if (H_num > 0) {
			final_name += "H" + H_num;
			
		}
		if (He_num > 0) {
			final_name += "He" + He_num;
			if (He_num > 0) {
				return false;
			}
		}
		if (Li_num > 0) {
			final_name += "Li" + Li_num;
			if (Li_num > 0) {
				return false;
			}
		}
		if (Be_num > 0) {
			final_name += "Be" + Be_num;
			if (Be_num > 0) {
				return false;
			}
		}
		if (B_num > 0) {
			final_name += "B" + B_num;
			if (B_num > 0) {
				return false;
			}
		}
		if (C_num > 0) {
			final_name += "C" + C_num;
			
		}
		if (N_num > 0) {
			final_name += "N" + N_num;
			
		}
		if (O_num > 0) {
			final_name += "O" + O_num;
			
		}
		if (F_num > 0) {
			final_name += "F" + F_num;
			
		}
		if (Ne_num > 0) {
			final_name += "Ne" + Ne_num;
			if (Ne_num > 0) {
				return false;
			}
		}
		if (Na_num > 0) {
			final_name += "Na" + Na_num;
			if (Na_num > 0) {
				return false;
			}
		}
		if (Mg_num > 0) {
			final_name += "Mg" + Mg_num;
			if (Mg_num > 0) {
				return false;
			}
		}
		if (Al_num > 0) {
			final_name += "Al" + Al_num;
			if (Al_num > 0) {
				return false;
			}
		}
		if (Si_num > 0) {
			final_name += "Si" + Si_num;
			if (Si_num > 0) {
				return false;
			}
		}
		if (P_num > 0) {
			final_name += "P" + P_num;
		}
		if (S_num > 0) {
			final_name += "S" + S_num;
		}
		if (Cl_num > 0) {
			final_name += "Cl" + Cl_num;
		}
		if (Ar_num > 0) {
			final_name += "Ar" + Ar_num;
			if (Ar_num > 0) {
				return false;
			}
		}
		if (K_num > 0) {
			final_name += "K" + K_num;
			if (K_num > 0) {
				return false;
			}
		}
		if (Ca_num > 0) {
			final_name += "Ca" + Ca_num;
			if (Ca_num > 0) {
				return false;
			}
		}
		if (Sc_num > 0) {
			final_name += "Sc" + Sc_num;
			if (Sc_num > 0) {
				return false;
			}
		}
		if (Ti_num > 0) {
			final_name += "Ti" + Ti_num;
			if (Ti_num > 0) {
				return false;
			}
		}
		if (V_num > 0) {
			final_name += "V" + V_num;
			if (V_num > 0) {
				return false;
			}
		}
		if (Cr_num > 0) {
			final_name += "Cr" + Cr_num;
			if (Cr_num > 0) {
				return false;
			}
		}
		if (Mn_num > 0) {
			final_name += "Mn" + Mn_num;
			if (Mn_num > 0) {
				return false;
			}
		}
		if (Fe_num > 0) {
			final_name += "Fe" + Fe_num;
			if (Fe_num > 0) {
				return false;
			}
		}
		if (Co_num > 0) {
			final_name += "Co" + Co_num;
			if (Co_num > 0) {
				return false;
			}
		}
		if (Ni_num > 0) {
			final_name += "Ni" + Ni_num;
			if (Ni_num > 0) {
				return false;
			}
		}
		if (Cu_num > 0) {
			final_name += "Cu" + Cu_num;
			if (Cu_num > 0) {
				return false;
			}
		}
		if (Zn_num > 0) {
			final_name += "Zn" + Zn_num;
			if (Zn_num > 0) {
				return false;
			}
		}
		if (Ga_num > 0) {
			final_name += "Ga" + Ga_num;
			if (Ga_num > 0) {
				return false;
			}
		}
		if (Ge_num > 0) {
			final_name += "Ge" + Ge_num;
			if (Ge_num > 0) {
				return false;
			}
		}
		if (As_num > 0) {
			final_name += "As" + As_num;
			if (As_num > 0) {
				return false;
			}
		}
		if (Se_num > 0) {
			final_name += "Se" + Se_num;
			if (Se_num > 0) {
				return false;
			}
		}
		if (Br_num > 0) {
			final_name += "Br" + Br_num;
		}
		if (Kr_num > 0) {
			final_name += "Kr" + Kr_num;
			if (Kr_num > 0) {
				return false;
			}
		}
		if (Rb_num > 0) {
			final_name += "Rb" + Rb_num;
			if (Rb_num > 0) {
				return false;
			}
		}
		if (Sr_num > 0) {
			final_name += "Sr" + Sr_num;
			if (Sr_num > 0) {
				return false;
			}
		}
		if (Y_num > 0) {
			final_name += "Y" + Y_num;
			if (Y_num > 0) {
				return false;
			}
		}
		if (Zr_num > 0) {
			final_name += "Zr" + Zr_num;
			if (Zr_num > 0) {
				return false;
			}
		}
		if (Nb_num > 0) {
			final_name += "Nb" + Nb_num;
			if (Nb_num > 0) {
				return false;
			}
		}
		if (Mo_num > 0) {
			final_name += "Mo" + Mo_num;
			if (Mo_num > 0) {
				return false;
			}
		}
		if (Tc_num > 0) {
			final_name += "Tc" + Tc_num;
			if (Tc_num > 0) {
				return false;
			}
		}
		if (Ru_num > 0) {
			final_name += "Ru" + Ru_num;
			if (Ru_num > 0) {
				return false;
			}
		}
		if (Rh_num > 0) {
			final_name += "Rh" + Rh_num;
			if (Rh_num > 0) {
				return false;
			}
		}
		if (Pd_num > 0) {
			final_name += "Pd" + Pd_num;
			if (Pd_num > 0) {
				return false;
			}
		}
		if (Ag_num > 0) {
			final_name += "Ag" + Ag_num;
			if (Ag_num > 0) {
				return false;
			}
		}
		if (Cd_num > 0) {
			final_name += "Cd" + Cd_num;
			if (Cd_num > 0) {
				return false;
			}
		}
		if (In_num > 0) {
			final_name += "In" + In_num;
			if (In_num > 0) {
				return false;
			}
		}
		if (Sn_num > 0) {
			final_name += "Sn" + Sn_num;
			if (Sn_num > 0) {
				return false;
			}
		}
		if (Sb_num > 0) {
			final_name += "Sb" + Sb_num;
			if (Sb_num > 0) {
				return false;
			}
		}
		if (Te_num > 0) {
			final_name += "Te" + Te_num;
			if (Te_num > 0) {
				return false;
			}
		}
		if (I_num > 0) {
			final_name += "I" + I_num;
			if (I_num > 0) {
				return false;
			}
		}
		if (Xe_num > 0) {
			final_name += "Xe" + Xe_num;
			if (Xe_num > 0) {
				return false;
			}
		}
		if (Cs_num > 0) {
			final_name += "Cs" + Cs_num;
			if (Cs_num > 0) {
				return false;
			}
		}
		if (Ba_num > 0) {
			final_name += "Ba" + Ba_num;
			if (Ba_num > 0) {
				return false;
			}
		}
		if (La_num > 0) {
			final_name += "La" + La_num;
			if (La_num > 0) {
				return false;
			}
		}
		if (Ce_num > 0) {
			final_name += "Ce" + Ce_num;
			if (Ce_num > 0) {
				return false;
			}
		}
		if (Pr_num > 0) {
			final_name += "Pr" + Pr_num;
			if (Pr_num > 0) {
				return false;
			}
		}
		if (Nd_num > 0) {
			final_name += "Nd" + Nd_num;
			if (Nd_num > 0) {
				return false;
			}
		}
		if (Pm_num > 0) {
			final_name += "Pm" + Pm_num;
			if (Pm_num > 0) {
				return false;
			}
		}
		if (Sm_num > 0) {
			final_name += "Sm" + Sm_num;
			if (Sm_num > 0) {
				return false;
			}
		}
		if (Eu_num > 0) {
			final_name += "Eu" + Eu_num;
			if (Eu_num > 0) {
				return false;
			}
		}
		if (Gd_num > 0) {
			final_name += "Gd" + Gd_num;
			if (Gd_num > 0) {
				return false;
			}
		}
		if (Tb_num > 0) {
			final_name += "Tb" + Tb_num;
			if (Tb_num > 0) {
				return false;
			}
		}
		if (Dy_num > 0) {
			final_name += "Dy" + Dy_num;
			if (Dy_num > 0) {
				return false;
			}
		}
		if (Ho_num > 0) {
			final_name += "Ho" + Ho_num;
			if (Ho_num > 0) {
				return false;
			}
		}
		if (Er_num > 0) {
			final_name += "Er" + Er_num;
			if (Er_num > 0) {
				return false;
			}
		}
		if (Tm_num > 0) {
			final_name += "Tm" + Tm_num;
			if (Tm_num > 0) {
				return false;
			}
		}
		if (Yb_num > 0) {
			final_name += "Yb" + Yb_num;
			if (Yb_num > 0) {
				return false;
			}
		}
		if (Lu_num > 0) {
			final_name += "Lu" + Lu_num;
			if (Lu_num > 0) {
				return false;
			}
		}
		if (Hf_num > 0) {
			final_name += "Hf" + Hf_num;
			if (Hf_num > 0) {
				return false;
			}
		}
		if (Ta_num > 0) {
			final_name += "Ta" + Ta_num;
			if (Ta_num > 0) {
				return false;
			}
		}
		if (W_num > 0) {
			final_name += "W" + W_num;
			if (W_num > 0) {
				return false;
			}
		}
		if (Re_num > 0) {
			final_name += "Re" + Re_num;
			if (Re_num > 0) {
				return false;
			}
		}
		if (Os_num > 0) {
			final_name += "Os" + Os_num;
			if (Os_num > 0) {
				return false;
			}
		}
		if (Ir_num > 0) {
			final_name += "Ir" + Ir_num;
			if (Ir_num > 0) {
				return false;
			}
		}
		if (Pt_num > 0) {
			final_name += "Pt" + Pt_num;
			if (Pt_num > 0) {
				return false;
			}
		}
		if (Au_num > 0) {
			final_name += "Au" + Au_num;
			if (Au_num > 0) {
				return false;
			}
		}
		if (Hg_num > 0) {
			final_name += "Hg" + Hg_num;
			if (Hg_num > 0) {
				return false;
			}
		}
		if (Tl_num > 0) {
			final_name += "Tl" + Tl_num;
			if (Tl_num > 0) {
				return false;
			}
		}
		if (Pb_num > 0) {
			final_name += "Pb" + Pb_num;
			if (Pb_num > 0) {
				return false;
			}
		}
		if (Bi_num > 0) {
			final_name += "Bi" + Bi_num;
			if (Bi_num > 0) {
				return false;
			}
		}
		if (Po_num > 0) {
			final_name += "Po" + Po_num;
			if (Po_num > 0) {
				return false;
			}
		}
		if (At_num > 0) {
			final_name += "At" + At_num;
			if (At_num > 0) {
				return false;
			}
		}
		if (Rn_num > 0) {
			final_name += "Rn" + Rn_num;
			if (Rn_num > 0) {
				return false;
			}
		}
		if (Fr_num > 0) {
			final_name += "Fr" + Fr_num;
			if (Fr_num > 0) {
				return false;
			}
		}
		if (Ra_num > 0) {
			final_name += "Ra" + Ra_num;
			if (Ra_num > 0) {
				return false;
			}
		}
		if (Ac_num > 0) {
			final_name += "Ac" + Ac_num;
			if (Ac_num > 0) {
				return false;
			}
		}
		if (Th_num > 0) {
			final_name += "Th" + Th_num;
			if (Th_num > 0) {
				return false;
			}
		}
		if (Pa_num > 0) {
			final_name += "Pa" + Pa_num;
			if (Pa_num > 0) {
				return false;
			}
		}
		if (U_num > 0) {
			final_name += "U" + U_num;
			if (U_num > 0) {
				return false;
			}
		}
		if (Np_num > 0) {
			final_name += "Np" + Np_num;
			if (Np_num > 0) {
				return false;
			}
		}
		if (Pu_num > 0) {
			final_name += "Pu" + Pu_num;
			if (Pu_num > 0) {
				return false;
			}
		}
		if (Am_num > 0) {
			final_name += "Am" + Am_num;
			if (Am_num > 0) {
				return false;
			}
		}
		if (Cm_num > 0) {
			final_name += "Cm" + Cm_num;
			if (Cm_num > 0) {
				return false;
			}
		}
		if (Bk_num > 0) {
			final_name += "Bk" + Bk_num;
			if (Bk_num > 0) {
				return false;
			}
		}
		if (Cf_num > 0) {
			final_name += "Cf" + Cf_num;
			if (Cf_num > 0) {
				return false;
			}
		}
		if (Es_num > 0) {
			final_name += "Es" + Es_num;
			if (Es_num > 0) {
				return false;
			}
		}
		if (Fm_num > 0) {
			final_name += "Fm" + Fm_num;
			if (Fm_num > 0) {
				return false;
			}
		}
		if (Md_num > 0) {
			final_name += "Md" + Md_num;
			if (Md_num > 0) {
				return false;
			}
		}
		if (No_num > 0) {
			final_name += "No" + No_num;
			if (No_num > 0) {
				return false;
			}
		}
		if (Lr_num > 0) {
			final_name += "Lr" + Lr_num;
			if (Lr_num > 0) {
				return false;
			}
		}
		if (Rf_num > 0) {
			final_name += "Rf" + Rf_num;
			if (Rf_num > 0) {
				return false;
			}
		}
		if (Db_num > 0) {
			final_name += "Db" + Db_num;
			if (Db_num > 0) {
				return false;
			}
		}
		if (Sg_num > 0) {
			final_name += "Sg" + Sg_num;
			if (Sg_num > 0) {
				return false;
			}
		}
		if (Bh_num > 0) {
			final_name += "Bh" + Bh_num;
			if (Bh_num > 0) {
				return false;
			}
		}
		if (Hs_num > 0) {
			final_name += "Hs" + Hs_num;
			if (Hs_num > 0) {
				return false;
			}
		}
		if (Mt_num > 0) {
			final_name += "Mt" + Mt_num;
			if (Mt_num > 0) {
				return false;
			}
		}
		if (Ds_num > 0) {
			final_name += "Ds" + Ds_num;
			if (Ds_num > 0) {
				return false;
			}
		}
		if (Rg_num > 0) {
			final_name += "Rg" + Rg_num;
			if (Rg_num > 0) {
				return false;
			}
		}
		if (Cn_num > 0) {
			final_name += "Cn" + Cn_num;
			if (Cn_num > 0) {
				return false;
			}
		}
		if (Uut_num > 0) {
			final_name += "Uut" + Uut_num;
			if (Uut_num > 0) {
				return false;
			}
		}
		if (Fl_num > 0) {
			final_name += "Fl" + Fl_num;
			if (Fl_num > 0) {
				return false;
			}
		}
		if (Uup_num > 0) {
			final_name += "Uup" + Uup_num;
			if (Uup_num > 0) {
				return false;
			}
		}
		if (Lv_num > 0) {
			final_name += "Lv" + Lv_num;
			if (Lv_num > 0) {
				return false;
			}
		}
		if (Uus_num > 0) {
			final_name += "Uus" + Uus_num;
			if (Uus_num > 0) {
				return false;
			}
		}
		if (Uuo_num > 0) {
			final_name += "Uuo" + Uuo_num;
			if (Uuo_num > 0) {
				return false;
			}
		}
		return true;
	}
	/**
	 * If the number of H, N, F, Cl, P, Br add up to even then satisfy the rule
	 * @return
	 */
	public static boolean check_hydrogen_rule(String formula) {
		
		int N_num = retrieve_num_element(formula, "N");
		int F_num = retrieve_num_element(formula, "F");
		int Cl_num = retrieve_num_element(formula, "Cl");
		int P_num = retrieve_num_element(formula, "P");
		int H_num = retrieve_num_element(formula, "H");
		int Br_num = retrieve_num_element(formula, "Br");
		int total = N_num + F_num + Cl_num + P_num + H_num + Br_num;
		if (total % 2 == 0) {
			return true;
		}
		return false;
	}
	/**
	 * This function only works for C, H, N, O, P, S, F, Cl, Br formulas
	 * @param formula
	 * @return
	 */
	public static String standardize_name(String formula) {
		if (!isNumeric(formula.substring(formula.length() - 1, formula.length()))) {
			formula = formula + "1";			
		}
		int Ci_num = retrieve_num_element(formula, "Ci");
		int H_num = retrieve_num_element(formula, "H");
		int He_num = retrieve_num_element(formula, "He");
		int Li_num = retrieve_num_element(formula, "Li");
		int Be_num = retrieve_num_element(formula, "Be");
		int B_num = retrieve_num_element(formula, "B");
		int C_num = retrieve_num_element(formula, "C");
		int N_num = retrieve_num_element(formula, "N");
		int O_num = retrieve_num_element(formula, "O");
		int F_num = retrieve_num_element(formula, "F");
		int Ne_num = retrieve_num_element(formula, "Ne");
		int Na_num = retrieve_num_element(formula, "Na");
		int Mg_num = retrieve_num_element(formula, "Mg");
		int Al_num = retrieve_num_element(formula, "Al");
		int Si_num = retrieve_num_element(formula, "Si");
		int P_num = retrieve_num_element(formula, "P");
		int S_num = retrieve_num_element(formula, "S");
		int Cl_num = retrieve_num_element(formula, "Cl");
		int Ar_num = retrieve_num_element(formula, "Ar");
		int K_num = retrieve_num_element(formula, "K");
		int Ca_num = retrieve_num_element(formula, "Ca");
		int Sc_num = retrieve_num_element(formula, "Sc");
		int Ti_num = retrieve_num_element(formula, "Ti");
		int V_num = retrieve_num_element(formula, "V");
		int Cr_num = retrieve_num_element(formula, "Cr");
		int Mn_num = retrieve_num_element(formula, "Mn");
		int Fe_num = retrieve_num_element(formula, "Fe");
		int Co_num = retrieve_num_element(formula, "Co");
		int Ni_num = retrieve_num_element(formula, "Ni");
		int Cu_num = retrieve_num_element(formula, "Cu");
		int Zn_num = retrieve_num_element(formula, "Zn");
		int Ga_num = retrieve_num_element(formula, "Ga");
		int Ge_num = retrieve_num_element(formula, "Ge");
		int As_num = retrieve_num_element(formula, "As");
		int Se_num = retrieve_num_element(formula, "Se");
		int Br_num = retrieve_num_element(formula, "Br");
		int Kr_num = retrieve_num_element(formula, "Kr");
		int Rb_num = retrieve_num_element(formula, "Rb");
		int Sr_num = retrieve_num_element(formula, "Sr");
		int Y_num = retrieve_num_element(formula, "Y");
		int Zr_num = retrieve_num_element(formula, "Zr");
		int Nb_num = retrieve_num_element(formula, "Nb");
		int Mo_num = retrieve_num_element(formula, "Mo");
		int Tc_num = retrieve_num_element(formula, "Tc");
		int Ru_num = retrieve_num_element(formula, "Ru");
		int Rh_num = retrieve_num_element(formula, "Rh");
		int Pd_num = retrieve_num_element(formula, "Pd");
		int Ag_num = retrieve_num_element(formula, "Ag");
		int Cd_num = retrieve_num_element(formula, "Cd");
		int In_num = retrieve_num_element(formula, "In");
		int Sn_num = retrieve_num_element(formula, "Sn");
		int Sb_num = retrieve_num_element(formula, "Sb");
		int Te_num = retrieve_num_element(formula, "Te");
		int I_num = retrieve_num_element(formula, "I");
		int Xe_num = retrieve_num_element(formula, "Xe");
		int Cs_num = retrieve_num_element(formula, "Cs");
		int Ba_num = retrieve_num_element(formula, "Ba");
		int La_num = retrieve_num_element(formula, "La");
		int Ce_num = retrieve_num_element(formula, "Ce");
		int Pr_num = retrieve_num_element(formula, "Pr");
		int Nd_num = retrieve_num_element(formula, "Nd");
		int Pm_num = retrieve_num_element(formula, "Pm");
		int Sm_num = retrieve_num_element(formula, "Sm");
		int Eu_num = retrieve_num_element(formula, "Eu");
		int Gd_num = retrieve_num_element(formula, "Gd");
		int Tb_num = retrieve_num_element(formula, "Tb");
		int Dy_num = retrieve_num_element(formula, "Dy");
		int Ho_num = retrieve_num_element(formula, "Ho");
		int Er_num = retrieve_num_element(formula, "Er");
		int Tm_num = retrieve_num_element(formula, "Tm");
		int Yb_num = retrieve_num_element(formula, "Yb");
		int Lu_num = retrieve_num_element(formula, "Lu");
		int Hf_num = retrieve_num_element(formula, "Hf");
		int Ta_num = retrieve_num_element(formula, "Ta");
		int W_num = retrieve_num_element(formula, "W");
		int Re_num = retrieve_num_element(formula, "Re");
		int Os_num = retrieve_num_element(formula, "Os");
		int Ir_num = retrieve_num_element(formula, "Ir");
		int Pt_num = retrieve_num_element(formula, "Pt");
		int Au_num = retrieve_num_element(formula, "Au");
		int Hg_num = retrieve_num_element(formula, "Hg");
		int Tl_num = retrieve_num_element(formula, "Tl");
		int Pb_num = retrieve_num_element(formula, "Pb");
		int Bi_num = retrieve_num_element(formula, "Bi");
		int Po_num = retrieve_num_element(formula, "Po");
		int At_num = retrieve_num_element(formula, "At");
		int Rn_num = retrieve_num_element(formula, "Rn");
		int Fr_num = retrieve_num_element(formula, "Fr");
		int Ra_num = retrieve_num_element(formula, "Ra");
		int Ac_num = retrieve_num_element(formula, "Ac");
		int Th_num = retrieve_num_element(formula, "Th");
		int Pa_num = retrieve_num_element(formula, "Pa");
		int U_num = retrieve_num_element(formula, "U");
		int Np_num = retrieve_num_element(formula, "Np");
		int Pu_num = retrieve_num_element(formula, "Pu");
		int Am_num = retrieve_num_element(formula, "Am");
		int Cm_num = retrieve_num_element(formula, "Cm");
		int Bk_num = retrieve_num_element(formula, "Bk");
		int Cf_num = retrieve_num_element(formula, "Cf");
		int Es_num = retrieve_num_element(formula, "Es");
		int Fm_num = retrieve_num_element(formula, "Fm");
		int Md_num = retrieve_num_element(formula, "Md");
		int No_num = retrieve_num_element(formula, "No");
		int Lr_num = retrieve_num_element(formula, "Lr");
		int Rf_num = retrieve_num_element(formula, "Rf");
		int Db_num = retrieve_num_element(formula, "Db");
		int Sg_num = retrieve_num_element(formula, "Sg");
		int Bh_num = retrieve_num_element(formula, "Bh");
		int Hs_num = retrieve_num_element(formula, "Hs");
		int Mt_num = retrieve_num_element(formula, "Mt");
		int Ds_num = retrieve_num_element(formula, "Ds");
		int Rg_num = retrieve_num_element(formula, "Rg");
		int Cn_num = retrieve_num_element(formula, "Cn");
		int Uut_num = retrieve_num_element(formula, "Uut");
		int Fl_num = retrieve_num_element(formula, "Fl");
		int Uup_num = retrieve_num_element(formula, "Uup");
		int Lv_num = retrieve_num_element(formula, "Lv");
		int Uus_num = retrieve_num_element(formula, "Uus");
		int Uuo_num = retrieve_num_element(formula, "Uuo");
		String final_name = "";
		if (Ci_num > 0) {
			final_name += "Ci" + Ci_num;
		}
		if (H_num > 0) {
			final_name += "H" + H_num;
		}
		if (He_num > 0) {
			final_name += "He" + He_num;
		}
		if (Li_num > 0) {
			final_name += "Li" + Li_num;
		}
		if (Be_num > 0) {
			final_name += "Be" + Be_num;
		}
		if (B_num > 0) {
			final_name += "B" + B_num;
		}
		if (C_num > 0) {
			final_name += "C" + C_num;
		}
		if (N_num > 0) {
			final_name += "N" + N_num;
		}
		if (O_num > 0) {
			final_name += "O" + O_num;
		}
		if (F_num > 0) {
			final_name += "F" + F_num;
		}
		if (Ne_num > 0) {
			final_name += "Ne" + Ne_num;
		}
		if (Na_num > 0) {
			final_name += "Na" + Na_num;
		}
		if (Mg_num > 0) {
			final_name += "Mg" + Mg_num;
		}
		if (Al_num > 0) {
			final_name += "Al" + Al_num;
		}
		if (Si_num > 0) {
			final_name += "Si" + Si_num;
		}
		if (P_num > 0) {
			final_name += "P" + P_num;
		}
		if (S_num > 0) {
			final_name += "S" + S_num;
		}
		if (Cl_num > 0) {
			final_name += "Cl" + Cl_num;
		}
		if (Ar_num > 0) {
			final_name += "Ar" + Ar_num;
		}
		if (K_num > 0) {
			final_name += "K" + K_num;
		}
		if (Ca_num > 0) {
			final_name += "Ca" + Ca_num;
		}
		if (Sc_num > 0) {
			final_name += "Sc" + Sc_num;
		}
		if (Ti_num > 0) {
			final_name += "Ti" + Ti_num;
		}
		if (V_num > 0) {
			final_name += "V" + V_num;
		}
		if (Cr_num > 0) {
			final_name += "Cr" + Cr_num;
		}
		if (Mn_num > 0) {
			final_name += "Mn" + Mn_num;
		}
		if (Fe_num > 0) {
			final_name += "Fe" + Fe_num;
		}
		if (Co_num > 0) {
			final_name += "Co" + Co_num;
		}
		if (Ni_num > 0) {
			final_name += "Ni" + Ni_num;
		}
		if (Cu_num > 0) {
			final_name += "Cu" + Cu_num;
		}
		if (Zn_num > 0) {
			final_name += "Zn" + Zn_num;
		}
		if (Ga_num > 0) {
			final_name += "Ga" + Ga_num;
		}
		if (Ge_num > 0) {
			final_name += "Ge" + Ge_num;
		}
		if (As_num > 0) {
			final_name += "As" + As_num;
		}
		if (Se_num > 0) {
			final_name += "Se" + Se_num;
		}
		if (Br_num > 0) {
			final_name += "Br" + Br_num;
		}
		if (Kr_num > 0) {
			final_name += "Kr" + Kr_num;
		}
		if (Rb_num > 0) {
			final_name += "Rb" + Rb_num;
		}
		if (Sr_num > 0) {
			final_name += "Sr" + Sr_num;
		}
		if (Y_num > 0) {
			final_name += "Y" + Y_num;
		}
		if (Zr_num > 0) {
			final_name += "Zr" + Zr_num;
		}
		if (Nb_num > 0) {
			final_name += "Nb" + Nb_num;
		}
		if (Mo_num > 0) {
			final_name += "Mo" + Mo_num;
		}
		if (Tc_num > 0) {
			final_name += "Tc" + Tc_num;
		}
		if (Ru_num > 0) {
			final_name += "Ru" + Ru_num;
		}
		if (Rh_num > 0) {
			final_name += "Rh" + Rh_num;
		}
		if (Pd_num > 0) {
			final_name += "Pd" + Pd_num;
		}
		if (Ag_num > 0) {
			final_name += "Ag" + Ag_num;
		}
		if (Cd_num > 0) {
			final_name += "Cd" + Cd_num;
		}
		if (In_num > 0) {
			final_name += "In" + In_num;
		}
		if (Sn_num > 0) {
			final_name += "Sn" + Sn_num;
		}
		if (Sb_num > 0) {
			final_name += "Sb" + Sb_num;
		}
		if (Te_num > 0) {
			final_name += "Te" + Te_num;
		}
		if (I_num > 0) {
			final_name += "I" + I_num;
		}
		if (Xe_num > 0) {
			final_name += "Xe" + Xe_num;
		}
		if (Cs_num > 0) {
			final_name += "Cs" + Cs_num;
		}
		if (Ba_num > 0) {
			final_name += "Ba" + Ba_num;
		}
		if (La_num > 0) {
			final_name += "La" + La_num;
		}
		if (Ce_num > 0) {
			final_name += "Ce" + Ce_num;
		}
		if (Pr_num > 0) {
			final_name += "Pr" + Pr_num;
		}
		if (Nd_num > 0) {
			final_name += "Nd" + Nd_num;
		}
		if (Pm_num > 0) {
			final_name += "Pm" + Pm_num;
		}
		if (Sm_num > 0) {
			final_name += "Sm" + Sm_num;
		}
		if (Eu_num > 0) {
			final_name += "Eu" + Eu_num;
		}
		if (Gd_num > 0) {
			final_name += "Gd" + Gd_num;
		}
		if (Tb_num > 0) {
			final_name += "Tb" + Tb_num;
		}
		if (Dy_num > 0) {
			final_name += "Dy" + Dy_num;
		}
		if (Ho_num > 0) {
			final_name += "Ho" + Ho_num;
		}
		if (Er_num > 0) {
			final_name += "Er" + Er_num;
		}
		if (Tm_num > 0) {
			final_name += "Tm" + Tm_num;
		}
		if (Yb_num > 0) {
			final_name += "Yb" + Yb_num;
		}
		if (Lu_num > 0) {
			final_name += "Lu" + Lu_num;
		}
		if (Hf_num > 0) {
			final_name += "Hf" + Hf_num;
		}
		if (Ta_num > 0) {
			final_name += "Ta" + Ta_num;
		}
		if (W_num > 0) {
			final_name += "W" + W_num;
		}
		if (Re_num > 0) {
			final_name += "Re" + Re_num;
		}
		if (Os_num > 0) {
			final_name += "Os" + Os_num;
		}
		if (Ir_num > 0) {
			final_name += "Ir" + Ir_num;
		}
		if (Pt_num > 0) {
			final_name += "Pt" + Pt_num;
		}
		if (Au_num > 0) {
			final_name += "Au" + Au_num;
		}
		if (Hg_num > 0) {
			final_name += "Hg" + Hg_num;
		}
		if (Tl_num > 0) {
			final_name += "Tl" + Tl_num;
		}
		if (Pb_num > 0) {
			final_name += "Pb" + Pb_num;
		}
		if (Bi_num > 0) {
			final_name += "Bi" + Bi_num;
		}
		if (Po_num > 0) {
			final_name += "Po" + Po_num;
		}
		if (At_num > 0) {
			final_name += "At" + At_num;
		}
		if (Rn_num > 0) {
			final_name += "Rn" + Rn_num;
		}
		if (Fr_num > 0) {
			final_name += "Fr" + Fr_num;
		}
		if (Ra_num > 0) {
			final_name += "Ra" + Ra_num;
		}
		if (Ac_num > 0) {
			final_name += "Ac" + Ac_num;
		}
		if (Th_num > 0) {
			final_name += "Th" + Th_num;
		}
		if (Pa_num > 0) {
			final_name += "Pa" + Pa_num;
		}
		if (U_num > 0) {
			final_name += "U" + U_num;
		}
		if (Np_num > 0) {
			final_name += "Np" + Np_num;
		}
		if (Pu_num > 0) {
			final_name += "Pu" + Pu_num;
		}
		if (Am_num > 0) {
			final_name += "Am" + Am_num;
		}
		if (Cm_num > 0) {
			final_name += "Cm" + Cm_num;
		}
		if (Bk_num > 0) {
			final_name += "Bk" + Bk_num;
		}
		if (Cf_num > 0) {
			final_name += "Cf" + Cf_num;
		}
		if (Es_num > 0) {
			final_name += "Es" + Es_num;
		}
		if (Fm_num > 0) {
			final_name += "Fm" + Fm_num;
		}
		if (Md_num > 0) {
			final_name += "Md" + Md_num;
		}
		if (No_num > 0) {
			final_name += "No" + No_num;
		}
		if (Lr_num > 0) {
			final_name += "Lr" + Lr_num;
		}
		if (Rf_num > 0) {
			final_name += "Rf" + Rf_num;
		}
		if (Db_num > 0) {
			final_name += "Db" + Db_num;
		}
		if (Sg_num > 0) {
			final_name += "Sg" + Sg_num;
		}
		if (Bh_num > 0) {
			final_name += "Bh" + Bh_num;
		}
		if (Hs_num > 0) {
			final_name += "Hs" + Hs_num;
		}
		if (Mt_num > 0) {
			final_name += "Mt" + Mt_num;
		}
		if (Ds_num > 0) {
			final_name += "Ds" + Ds_num;
		}
		if (Rg_num > 0) {
			final_name += "Rg" + Rg_num;
		}
		if (Cn_num > 0) {
			final_name += "Cn" + Cn_num;
		}
		if (Uut_num > 0) {
			final_name += "Uut" + Uut_num;
		}
		if (Fl_num > 0) {
			final_name += "Fl" + Fl_num;
		}
		if (Uup_num > 0) {
			final_name += "Uup" + Uup_num;
		}
		if (Lv_num > 0) {
			final_name += "Lv" + Lv_num;
		}
		if (Uus_num > 0) {
			final_name += "Uus" + Uus_num;
		}
		if (Uuo_num > 0) {
			final_name += "Uuo" + Uuo_num;
		}
		return final_name;
	}
	/**
	 * Parsing the string and retrieve the number of elements
	 * @param str
	 * @param query_element
	 * @return
	 */
	public static int retrieve_num_element(String str, String query_element) {
		int total = 0;
		for (int i = 0; i < str.length() - (query_element.length() - 1); i++) {
			if (str.substring(i, i + query_element.length()).equals(query_element)) {
				boolean verify = true;
				if (query_element.length() == 1 && i + 1 < str.length()) {
					if (!(str.substring(i + 1, i + 2).matches("[0-9]") || Character.isUpperCase(str.substring(i + 1, i + 2).charAt(0))   )) {
						verify = false;
					}
				}
				if (verify) {
					int value = 0;
					String num_str = "";
					for (int j = i + query_element.length(); j < str.length(); j++) {
						
						if (str.substring(j, j + 1).matches("[0-9]")) {
							num_str += str.substring(j, j + 1);
							value = new Integer(num_str);
						} else {
							
							if (value == 0) {
								value = 1;
							}
							break;
						}
					}
					if (i == str.length() - 1) {
						value = 1;
					}
					total += value;
				}
			}
		}
		return total;
	}
	public static double getMonoisotopicMass(String formula) {
		if (!isNumeric(formula.substring(formula.length() - 1, formula.length()))) {
			formula = formula + "1";			
		}
		int Ci_num = retrieve_num_element(formula, "Ci");
		int H_num = retrieve_num_element(formula, "H");
		int He_num = retrieve_num_element(formula, "He");
		int Li_num = retrieve_num_element(formula, "Li");
		int Be_num = retrieve_num_element(formula, "Be");
		int B_num = retrieve_num_element(formula, "B");
		int C_num = retrieve_num_element(formula, "C");
		int N_num = retrieve_num_element(formula, "N");
		int O_num = retrieve_num_element(formula, "O");
		int F_num = retrieve_num_element(formula, "F");
		int Ne_num = retrieve_num_element(formula, "Ne");
		int Na_num = retrieve_num_element(formula, "Na");
		int Mg_num = retrieve_num_element(formula, "Mg");
		int Al_num = retrieve_num_element(formula, "Al");
		int Si_num = retrieve_num_element(formula, "Si");
		int P_num = retrieve_num_element(formula, "P");
		int S_num = retrieve_num_element(formula, "S");
		int Cl_num = retrieve_num_element(formula, "Cl");
		int Ar_num = retrieve_num_element(formula, "Ar");
		int K_num = retrieve_num_element(formula, "K");
		int Ca_num = retrieve_num_element(formula, "Ca");
		int Sc_num = retrieve_num_element(formula, "Sc");
		int Ti_num = retrieve_num_element(formula, "Ti");
		int V_num = retrieve_num_element(formula, "V");
		int Cr_num = retrieve_num_element(formula, "Cr");
		int Mn_num = retrieve_num_element(formula, "Mn");
		int Fe_num = retrieve_num_element(formula, "Fe");
		int Co_num = retrieve_num_element(formula, "Co");
		int Ni_num = retrieve_num_element(formula, "Ni");
		int Cu_num = retrieve_num_element(formula, "Cu");
		int Zn_num = retrieve_num_element(formula, "Zn");
		int Ga_num = retrieve_num_element(formula, "Ga");
		int Ge_num = retrieve_num_element(formula, "Ge");
		int As_num = retrieve_num_element(formula, "As");
		int Se_num = retrieve_num_element(formula, "Se");
		int Br_num = retrieve_num_element(formula, "Br");
		int Kr_num = retrieve_num_element(formula, "Kr");
		int Rb_num = retrieve_num_element(formula, "Rb");
		int Sr_num = retrieve_num_element(formula, "Sr");
		int Y_num = retrieve_num_element(formula, "Y");
		int Zr_num = retrieve_num_element(formula, "Zr");
		int Nb_num = retrieve_num_element(formula, "Nb");
		int Mo_num = retrieve_num_element(formula, "Mo");
		int Tc_num = retrieve_num_element(formula, "Tc");
		int Ru_num = retrieve_num_element(formula, "Ru");
		int Rh_num = retrieve_num_element(formula, "Rh");
		int Pd_num = retrieve_num_element(formula, "Pd");
		int Ag_num = retrieve_num_element(formula, "Ag");
		int Cd_num = retrieve_num_element(formula, "Cd");
		int In_num = retrieve_num_element(formula, "In");
		int Sn_num = retrieve_num_element(formula, "Sn");
		int Sb_num = retrieve_num_element(formula, "Sb");
		int Te_num = retrieve_num_element(formula, "Te");
		int I_num = retrieve_num_element(formula, "I");
		int Xe_num = retrieve_num_element(formula, "Xe");
		int Cs_num = retrieve_num_element(formula, "Cs");
		int Ba_num = retrieve_num_element(formula, "Ba");
		int La_num = retrieve_num_element(formula, "La");
		int Ce_num = retrieve_num_element(formula, "Ce");
		int Pr_num = retrieve_num_element(formula, "Pr");
		int Nd_num = retrieve_num_element(formula, "Nd");
		int Pm_num = retrieve_num_element(formula, "Pm");
		int Sm_num = retrieve_num_element(formula, "Sm");
		int Eu_num = retrieve_num_element(formula, "Eu");
		int Gd_num = retrieve_num_element(formula, "Gd");
		int Tb_num = retrieve_num_element(formula, "Tb");
		int Dy_num = retrieve_num_element(formula, "Dy");
		int Ho_num = retrieve_num_element(formula, "Ho");
		int Er_num = retrieve_num_element(formula, "Er");
		int Tm_num = retrieve_num_element(formula, "Tm");
		int Yb_num = retrieve_num_element(formula, "Yb");
		int Lu_num = retrieve_num_element(formula, "Lu");
		int Hf_num = retrieve_num_element(formula, "Hf");
		int Ta_num = retrieve_num_element(formula, "Ta");
		int W_num = retrieve_num_element(formula, "W");
		int Re_num = retrieve_num_element(formula, "Re");
		int Os_num = retrieve_num_element(formula, "Os");
		int Ir_num = retrieve_num_element(formula, "Ir");
		int Pt_num = retrieve_num_element(formula, "Pt");
		int Au_num = retrieve_num_element(formula, "Au");
		int Hg_num = retrieve_num_element(formula, "Hg");
		int Tl_num = retrieve_num_element(formula, "Tl");
		int Pb_num = retrieve_num_element(formula, "Pb");
		int Bi_num = retrieve_num_element(formula, "Bi");
		int Po_num = retrieve_num_element(formula, "Po");
		int At_num = retrieve_num_element(formula, "At");
		int Rn_num = retrieve_num_element(formula, "Rn");
		int Fr_num = retrieve_num_element(formula, "Fr");
		int Ra_num = retrieve_num_element(formula, "Ra");
		int Ac_num = retrieve_num_element(formula, "Ac");
		int Th_num = retrieve_num_element(formula, "Th");
		int Pa_num = retrieve_num_element(formula, "Pa");
		int U_num = retrieve_num_element(formula, "U");
		int Np_num = retrieve_num_element(formula, "Np");
		int Pu_num = retrieve_num_element(formula, "Pu");
		int Am_num = retrieve_num_element(formula, "Am");
		int Cm_num = retrieve_num_element(formula, "Cm");
		int Bk_num = retrieve_num_element(formula, "Bk");
		int Cf_num = retrieve_num_element(formula, "Cf");
		int Es_num = retrieve_num_element(formula, "Es");
		int Fm_num = retrieve_num_element(formula, "Fm");
		int Md_num = retrieve_num_element(formula, "Md");
		int No_num = retrieve_num_element(formula, "No");
		int Lr_num = retrieve_num_element(formula, "Lr");
		int Rf_num = retrieve_num_element(formula, "Rf");
		int Db_num = retrieve_num_element(formula, "Db");
		int Sg_num = retrieve_num_element(formula, "Sg");
		int Bh_num = retrieve_num_element(formula, "Bh");
		int Hs_num = retrieve_num_element(formula, "Hs");
		int Mt_num = retrieve_num_element(formula, "Mt");
		int Ds_num = retrieve_num_element(formula, "Ds");
		int Rg_num = retrieve_num_element(formula, "Rg");
		int Cn_num = retrieve_num_element(formula, "Cn");
		int Uut_num = retrieve_num_element(formula, "Uut");
		int Fl_num = retrieve_num_element(formula, "Fl");
		int Uup_num = retrieve_num_element(formula, "Uup");
		int Lv_num = retrieve_num_element(formula, "Lv");
		int Uus_num = retrieve_num_element(formula, "Uus");
		int Uuo_num = retrieve_num_element(formula, "Uuo");
		double C_MASS = 12;
		double Ci_MASS = 13.0033548378;
		double H_MASS = 1.0078250321;
		double O_MASS = 15.9949146221;
		double P_MASS = 30.97376151;
		double N_MASS = 14.0030740052;
		double S_MASS = 31.97207069;
		double Cl_MASS = 34.96885272;
		double Br_MASS = 78.9183361;
		double F_MASS = 18.99840322;
		return Ci_num * Ci_MASS + H_num * H_MASS + C_num * C_MASS + O_num * O_MASS 
				+ N_num * N_MASS + P_num * P_MASS + S_num * S_MASS
				+ Cl_num * Cl_MASS + Br_num * Br_MASS + F_num * F_MASS;
		
	}
}
