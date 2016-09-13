package SIMULATE_FORMULA;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.formula.IsotopeContainer;
import org.openscience.cdk.formula.IsotopePattern;
import org.openscience.cdk.formula.IsotopePatternGenerator;
import org.openscience.cdk.formula.IsotopePatternSimilarity;
import org.openscience.cdk.formula.MassToFormulaTool;
import org.openscience.cdk.formula.MolecularFormulaRange;
import org.openscience.cdk.formula.rules.ChargeRule;
import org.openscience.cdk.formula.rules.ElementRule;
import org.openscience.cdk.formula.rules.IsotopePatternRule;
import org.openscience.cdk.formula.rules.MMElementRule;
import org.openscience.cdk.formula.rules.NitrogenRule;
import org.openscience.cdk.formula.rules.RDBERule;
import org.openscience.cdk.formula.rules.ToleranceRangeRule;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecularFormulaSet;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.Atom;
import org.openscience.cdk.Isotope;
import org.openscience.cdk.formula.MolecularFormula;
import org.openscience.cdk.DefaultChemObjectBuilder;



public class SimulateFormula {
	public static int C_MIN = 0;
	public static int C_MAX = 105;
	public static int N_MIN = 0;
	public static int N_MAX = 30;
	public static int O_MIN = 0;
	public static int O_MAX = 40;
	public static int H_MIN = 0;
	public static int H_MAX = 170;
	public static int S_MIN = 0;
	public static int S_MAX = 4;
	public static int P_MIN = 0;
	public static int P_MAX = 4;
        public static int Cl_MIN = 0;
        public static int Cl_MAX = 4;
        public static int Br_MIN = 0;
        public static int Br_MAX = 4;
        public static int F_MIN = 0;
        public static int F_MAX = 4;

	public static double _tolerance = 0.5; //ppm
	public static String _SP = "";
	public static void main(String[] args) {
		try {
			_tolerance = new Double(args[1]);
			
			
			HashMap formula_db = new HashMap();
			String fileName = args[0];
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			in.readLine();
			while (in.ready()) {
				
				String str = in.readLine();
				String[] split = str.split("\t");
				
				String compound_mass = split[0];
				if (!compound_mass.equals("")) {
					//masses_hmdb.add(new Double(split[3]));
					formula_db.put(compound_mass, compound_mass);
				}
			}
			in.close();
			
			HashMap sampled_mass = new HashMap();
			//FileWriter fwriter = new FileWriter("C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Formula_Expansion");
			FileWriter fwriter = new FileWriter(args[2]); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Formula_Expansion");
			BufferedWriter out = new BufferedWriter(fwriter);
			
			//FileWriter fwriter2 = new FileWriter("C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\All_Formulas.txt");
			FileWriter fwriter2 = new FileWriter(args[3]); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\All_Formulas.txt");
			BufferedWriter out2 = new BufferedWriter(fwriter2);
			
			String value = args[4];
			_SP = value;
			/*if (value.equals("true")) {
				_SP = true;
			} else {
				_SP = false;
			}*/
			
			//int i = 100;
			//for (int i = 0; i < 1000; i++) {
			
			
			Iterator itr = formula_db.keySet().iterator();
			while (itr.hasNext()) {
				double reference_compound_mass = new Double((String)itr.next());
				double mass_db = reference_compound_mass;
				if (mass_db < 1501) {
					if (!sampled_mass.containsKey(mass_db)) {
						IMolecularFormulaSet mol_set = getFormulas(mass_db);
						if (mol_set != null) {
							
							System.out.println(mass_db + "\t" + mol_set.size());
							//out.write(mass_hmdb + "\t" + mol_set.size() + "\n");
							//out.flush();
							
							//out2.write(mass_hmdb + "");
							
							String[] formulas = getFormulaMasses(mol_set);
							sampled_mass.put(mass_db, formulas);
							//for (String formula: formulas) {
							//	out2.write("\t" + formula);
							//}
							//out2.write("\n");
							//out2.flush();
						}
					}
					if (sampled_mass.containsKey(mass_db)) {
						String[] formulas = (String[])sampled_mass.get(mass_db);
						out.write(mass_db + "\t" + formulas.length + "\n");
						out.flush();
						
						out2.write(reference_compound_mass + "\t" + mass_db + "");
						for (String formula: formulas) {
							out2.write("\t" + formula);
						}
						out2.write("\n");
						out2.flush();	
					}
					
				}
				//
			}
			out.close();
			out2.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public static String[] getFormulas(IMolecularFormulaSet mfSet) {
		String[] formulas = new String[mfSet.size()];
		int i = 0;
		Iterable<IMolecularFormula> iterator = mfSet.molecularFormulas();
		Iterator itr = iterator.iterator(); 
		while (itr.hasNext()) {
			MolecularFormula formula = (MolecularFormula)itr.next();						
			formulas[i] = getFormula(formula);
			i++;
		}
		return formulas;
	}
	public static String[] getFormulaMasses(IMolecularFormulaSet mfSet) {
		String[] formulas = new String[mfSet.size()];
		int i = 0;
		Iterable<IMolecularFormula> iterator = mfSet.molecularFormulas();
		Iterator itr = iterator.iterator(); 
		while (itr.hasNext()) {
			MolecularFormula formula = (MolecularFormula)itr.next();						
			formulas[i] = getFormulaMass(formula);
			i++;
		}
		return formulas;
	}
	public static String getFormula(IMolecularFormula formula) {
		
		String formulastr = "";
		
		Iterable<IIsotope> iterator2 = formula.isotopes();
		Iterator itr2 = iterator2.iterator();
		while (itr2.hasNext()) {
			Isotope isotope = (Isotope)itr2.next();
			
			formulastr += isotope.getSymbol() + formula.getIsotopeCount(isotope);
			//System.out.print(isotope.getSymbol() + formula.getIsotopeCount(isotope));
			//isotope.exactMass + "\t"
							
		}
		return formulastr;
	}
public static String getFormulaMass(IMolecularFormula formula) {
		
		String formulastr = "";
		IAtomContainer atomContainer = MolecularFormulaManipulator.getAtomContainer(formula);
		//double mass = AtomContainerManipulator.getNaturalExactMass(atomContainer);
		double mass = AtomContainerManipulator.getTotalExactMass(atomContainer);
		Iterable<IIsotope> iterator2 = formula.isotopes();
		Iterator itr2 = iterator2.iterator();
		while (itr2.hasNext()) {
			Isotope isotope = (Isotope)itr2.next();
			
			formulastr += isotope.getSymbol() + formula.getIsotopeCount(isotope);
			//System.out.print(isotope.getSymbol() + formula.getIsotopeCount(isotope));
			//isotope.exactMass + "\t"
							
		}
		return formulastr + ":" + mass;
	}
	public static IMolecularFormulaSet getFormulas(double myMass) throws CDKException {
		MassToFormulaTool mf = new MassToFormulaTool(DefaultChemObjectBuilder.getInstance());		
		// tolerance rule
		
		//double diff = _tolerance * myMass / 1E6;
		double diff = _tolerance;
		
		ToleranceRangeRule tolerance = new ToleranceRangeRule();
		Object[] objs = new Object[2];
		objs[0] = new Double(0.0);
		objs[1] = new Double(diff);
		tolerance.setParameters(objs);
		
		// nitrogen rule
		NitrogenRule nitrogen = new NitrogenRule();
		
		//RDBRule
		RDBERule rdbrule = new RDBERule();
		
		//Element Rule
		ElementRule elementRule = new ElementRule();		
		//IsotopePattern isotopes = new IsotopePattern();				
		
		MolecularFormulaRange atoms = new MolecularFormulaRange();
		Atom C = new Atom("C");
		Atom N = new Atom("N");
		Atom S = new Atom("S");
		Atom O = new Atom("O");
		Atom H = new Atom("H");
		Atom P = new Atom("P");
		Atom Br = new Atom("Br");
		Atom Cl = new Atom("Cl");
		Atom F = new Atom("F");
		
		IsotopeFactory f;
		try {

		
			f = IsotopeFactory.getInstance(C.getBuilder());
			f.configure(C);
			atoms.addIsotope(C, C_MIN, C_MAX);
			
			f = IsotopeFactory.getInstance(N.getBuilder());
			f.configure(N);			
			atoms.addIsotope(N, N_MIN, N_MAX);
			
			if (_SP.equals("S")) {
				f = IsotopeFactory.getInstance(S.getBuilder());
				f.configure(S);			
				atoms.addIsotope(S, S_MIN, S_MAX);
			} else if (_SP.equals("P")) {	
				f = IsotopeFactory.getInstance(P.getBuilder());
				f.configure(P);			
				atoms.addIsotope(P, P_MIN, P_MAX);
				
			} else if (_SP.equals("F")) {
                                f = IsotopeFactory.getInstance(F.getBuilder());
                                f.configure(F);
                                atoms.addIsotope(F, F_MIN, F_MAX);				
			} else if (_SP.equals("Br")) {
                                f = IsotopeFactory.getInstance(Br.getBuilder());
                                f.configure(Br);
                                atoms.addIsotope(Br, Br_MIN, Br_MAX);
			} else if (_SP.equals("Cl")) {
                                f = IsotopeFactory.getInstance(Cl.getBuilder());
                                f.configure(Cl);
                                atoms.addIsotope(Cl, Cl_MIN, Cl_MAX);
			}
			
			f = IsotopeFactory.getInstance(O.getBuilder());
			f.configure(O);			
			atoms.addIsotope(O, O_MIN, O_MAX);
			
			f = IsotopeFactory.getInstance(H.getBuilder());
			f.configure(H);			
			atoms.addIsotope(H, H_MIN, H_MAX);

			//System.out.println(C_MIN + "\t" + C_MAX);
			//System.out.println(O_MIN + "\t" + O_MAX);
			//System.out.println(N_MIN + "\t" + N_MAX);
			//System.out.println(S_MIN + "\t" + S_MAX);
			//System.out.println(H_MIN + "\t" + H_MAX);
			
			/*
			f = IsotopeFactory.getInstance(C.getBuilder());
			f.configure(C);
			atoms.addIsotope(C, 0, 30);
			
			f = IsotopeFactory.getInstance(N.getBuilder());
			f.configure(N);			
			atoms.addIsotope(N, 0, 10);
			
			f = IsotopeFactory.getInstance(S.getBuilder());
			f.configure(S);			
			atoms.addIsotope(S, 0, 10);
			
			f = IsotopeFactory.getInstance(O.getBuilder());
			f.configure(O);			
			atoms.addIsotope(O, 0, 15);
			
			f = IsotopeFactory.getInstance(H.getBuilder());
			f.configure(H);			
			atoms.addIsotope(H, 0, 60);
			*/
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}					
		
	 	objs = new Object[1];
		objs[0] = atoms;
		elementRule.setParameters(objs);
		
		
		MMElementRule mmelement = new MMElementRule();
		objs = new Object[2];
		//objs[0] = MMElementRule.Database.WILEY;
		objs[0] = MMElementRule.Database.DNP;
		objs[1] = MMElementRule.RangeMass.Minus2000;
		mmelement.setParameters(objs);
		
		
		
		IsotopePatternRule isotopePattern = new IsotopePatternRule(); // 
		//objs = new Object[2];
		//List<double[]> isotopes = testIsotopeContainer();
		//objs[0] = isotopes;		
		//objs[1] = new Double(0.01);
		//isotopePattern.setParameters(objs);
		ChargeRule chargeRule = new ChargeRule(); // default is 0 neutral charge
		
		LinkedList list = new LinkedList();
		list.add(mmelement);
		list.add(tolerance);
		list.add(nitrogen);
		list.add(elementRule);
		list.add(rdbrule);
		list.add(isotopePattern);
		list.add(chargeRule);
		
		mf.setRestrictions(list);		
		//mf.setDefaultRestrictions();		
		//mf.setRestrictions(arg0);;
						
		IMolecularFormulaSet mfSet = mf.generate(myMass);
		if (mfSet == null) {
			return null;
		} else {
			
			//printFormulas(mfSet);
		}
		return mfSet;
	}
	public static void printFormulas(IMolecularFormulaSet mfSet) {
		Iterable<IMolecularFormula> iterator = mfSet.molecularFormulas();
		Iterator itr = iterator.iterator(); 
		while (itr.hasNext()) {
			MolecularFormula formula = (MolecularFormula)itr.next();						
			
			Iterable<IIsotope> iterator2 = formula.isotopes();
			Iterator itr2 = iterator2.iterator();
			while (itr2.hasNext()) {
				Isotope isotope = (Isotope)itr2.next();
				System.out.print(isotope.getSymbol() + formula.getIsotopeCount(isotope));
				//isotope.exactMass + "\t"
								
			}
			
			//System.out.println();
			
			IsotopePatternGenerator isotopeGenerator = new IsotopePatternGenerator(0.1);
			IsotopePattern pattern = isotopeGenerator.getIsotopes(formula);
			IsotopePattern reference_pattern = testIsotopeContainer();
			IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
			
			
			System.out.println("Num isotope: " + pattern.getNumberOfIsotopes());
			//printIsotope(pattern);
			//for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {
				//IsotopeContainer container = pattern.getIsotope(i);
				//System.out.println(container.getMass() + "\t" + container.getIntensity());				
			//}
			System.out.println("Simiarlity: " + similarity.compare(pattern, reference_pattern));
			//System.out.println(formula);	
			
		}
		System.out.println(mfSet.size());
	}
	//public static LinkedList<double[]> testIsotopeContainer() {
		public static IsotopePattern testIsotopeContainer() {
			/* 
			301 100
			302 7.8
			303 7.1
			304 0.4
			305 0.1
			*/
			//C6H5O12S1
			LinkedList<double[]> list = new LinkedList();
			double[] d = new double[2];
			d[0] = 300.9501716;;
			d[1] = 1.0;
			list.add(d);
			d[0] = 301.95352643999996;
			d[1] = 0.06489436975639341;
			list.add(d);
			d[0] = 302.9459675;
			d[1] = 0.0451911935110081;		
			list.add(d);
			d[0] = 302.95441798;
			d[1] = 0.024659923614382958;		
			list.add(d);

			IsotopePattern result = new IsotopePattern();
			result.addIsotope(new IsotopeContainer(300.9501716, 1.0));
			result.addIsotope(new IsotopeContainer(301.95352643999996, 0.06489436975639341));
			result.addIsotope(new IsotopeContainer(302.9459675, 0.0451911935110081));
			result.addIsotope(new IsotopeContainer(302.95441798, 0.024659923614382958));
			
			//list.add(result);
			return result;
		}
}

