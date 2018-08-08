/************************************************************           
 * MutPanning - Step 10										*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2018 									*   
 *															*   
 * License:		Public Domain								*   
 *															*   
 * Summary: This script optimizes the parameters of the 	*
 * CBASE model to model the background distribution of		*
 * synonymous mutations across the exome. For this purpose,	*
 * we use a quasi-Newton approach (Limited-Memory BFGS) for	*
 * the optimization. For each model we start with 150 random*   
 * initial choices for the parameters, perform L-BFGS and 	*
 * select the best solution.								*
 *************************************************************/

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.apache.commons.math3.special.Gamma;
import jdistlib.math.Bessel;

public class CBASE_Solutions {
	static String[] index_header_samples={"ID","Sample","Cohort"};
	
	/*
	 * argument0 root file
	 * arugment1 sample file
	 */
	
	public static void main (String[] args){
		//Determine all entity names, as clustering is performed separately for each Step this is needed to coordinate the order
		String[] entity_name=new String[0];
		try{
			FileInputStream in=new FileInputStream(args[1]);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			String s="";
			ArrayList<String> aa =new ArrayList<String>();
			while((s=input.readLine())!=null){
				String e=s.split("	")[index_header[2]];
				if(!contains(e,aa)){
					aa.add(e);
				}
			}
			input.close();
			Collections.sort(aa);
			aa.add("PanCancer");
			entity_name=new String[aa.size()];
			for (int i=0;i<aa.size();i++){
				entity_name[i]=aa.get(i);
			}
		}
		catch(Exception e){
			System.out.println(e);
		}
		
		Comparator<int[]> comp=(int[] a, int[] b)->{
			return new Integer(a[0]).compareTo(b[0]);
		};
		
		if(!new File(args[0]+"CBASE/Parameters_Summary/").exists()){
			new File(args[0]+"CBASE/Parameters_Summary/").mkdir();
		}
		
		for (int k=0;k<entity_name.length;k++){
			System.out.println(entity_name[k]);
			ArrayList<Integer> counts=new ArrayList<Integer>();
			try{
				FileInputStream in=new FileInputStream(args[0]+"CBASE/Counts/CountSilent"+entity_name[k]+".txt");
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				input.readLine();
				String s="";
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					if(!is_olfactory(t[0])){
						counts.add(Integer.parseInt(t[2]));
					}
				}
				input.close();
				
				
			}
			catch(Exception e){
				System.out.println(e);
			}
			ArrayList<int[]> histogram=histogram(counts);
			Collections.sort(histogram,comp);
			double[][] solution=minimize_parallel(150, histogram,  24, 100);
			
			try{
				FileWriter out=new FileWriter(args[0]+"CBASE/Parameters_Summary/Parameters"+entity_name[k]+".txt");
				BufferedWriter output=new BufferedWriter(out);
				for (int i=0;i<solution.length;i++){
					output.write(ou(solution[i])+", "+(i+1));
					output.newLine();
				}
				output.close();
			}
			catch(Exception e){
				System.out.println(e);
			}
		}
	
	}
	
	public static int[] index_header(String[] header, String[] ideal_header){
		int[] indices=new int[ideal_header.length];
		for (int i=0;i<ideal_header.length;i++){
			int index=-1;
			for (int j=0;j<header.length;j++){
				if(header[j].equals(ideal_header[i])){
					index=j;
					break;
				}
			}
			indices[i]=index;
		}
		return indices;
	}
	
	public static boolean contains (String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
	
	public static double[] improve(int mod_C,  ArrayList<int[]> histogram, double[] start, double[][] bounds, int iterations){
		double[] a=minimize_neg_ln_L(discard_last(start),mod_C,bounds,histogram,iterations);
		double value=model(histogram,discard_last(a),mod_C);
		for (int i=0;i<100;i++){
			a=minimize_neg_ln_L(discard_last(a),mod_C,bounds,histogram,iterations);
			double value_new=model(histogram,discard_last(a),mod_C);
			if(value==value_new){
				break;
			}
			value=value_new;
			//System.out.println(ou(a));
		}
		return a;
	}
	public static double[] discard_last(double[] a){
		double[] b=new double[a.length-1];
		for (int i=0;i<a.length-1;i++){
			b[i]=a[i];
		}
		return b;
	}
	
	public static double [] minimize_neg_ln_L(double[] start, int mod, double[][] bounds, ArrayList<int[]> histo ,int iterations){//
		//bounds are currently ignored
		double[][][] standard_bounds={{{0,100},{0,100}},{{0,100},{0,100}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}}};
		
		
		//double[] x=start;
		//double[] g=gradient(x);
		ArrayList<double[]> g=new ArrayList<double[]>();//stores the last m+1 entries
		ArrayList<double[]> x=new ArrayList<double[]>();//stores the last m+1 entries
		
		x.add(start);
		g.add(gradient(x.get(x.size()-1),mod, histo));
//		System.out.println((0)+"	"+ou(x.get(x.size()-1))+"	"+model(histo,x.get(x.size()-1),mod));
		
		double value_best_solution=model(histo,start,mod);
		double[] best_solution=clone(start);
		
		outer:
		for (int iter=0;iter<iterations;iter++){
			double[] z=new double[0];
			//double scaling=1;
			if(x.size()>1){
				double[] q=clone(g.get(g.size()-1));
				double [] alpha=new double[x.size()-1];
				for (int i=1;i<x.size();i++){
					double[] y=diff(g.get(x.size()-i),g.get(x.size()-i-1));
					double[] s=diff(x.get(x.size()-i),x.get(x.size()-i-1));
					double rho=1/product_tv_v(y,s);
					alpha[x.size()-i-1]=rho*product_tv_v(s,q);
					q=diff(q,product(alpha[i-1],y));
				}
				
				
				double[] y=diff(g.get(g.size()-1),g.get(g.size()-2));
				double[] s=diff(x.get(x.size()-1),x.get(x.size()-2));
				
				//double scaling=product_tv_v(s,y)/sq_norm(y);
				double[][] h_k_0=product(1.0/sq_norm(y),product_v_tv(y,s));
				z=product(h_k_0,q);
				for (int i=1;i<x.size();i++){
					y=diff(g.get(i),g.get(i-1));
					s=diff(x.get(i),x.get(i-1));
					double rho=1/product_tv_v(y,s);
					double beta=rho*product_tv_v(y,z);
					z=add(z,product(alpha[i-1]-beta,s));
				}
				//System.out.println("z: "+ou(z));
				if(is_nan(z)){
					break outer;
				}
				
				z=product(1/Math.sqrt(sq_norm(z)),z);
				
				double min=Double.MAX_VALUE;
				double k_min=model(histo,x.get(x.size()-1),mod);
				for (double k=0.01;k<=1;k+=0.01){
					double m=model(histo,diff(x.get(x.size()-1),product(k,z)),mod);
					if(m<min&&in_bound(diff(x.get(x.size()-1),product(k,z)),bounds)){
						min=m;
						k_min=k;
					}
				}
				if(min!=Double.MAX_VALUE){
					z=product(k_min,z);
				}
				else{
					min=Double.MAX_VALUE;
					k_min=model(histo,x.get(x.size()-1),mod);
					for (double k=0.0001;k<=0.01;k+=0.0001){
						double m=model(histo,diff(x.get(x.size()-1),product(k,z)),mod);
						if(m<min&&in_bound(diff(x.get(x.size()-1),product(k,z)),bounds)){
							min=m;
							k_min=k;
						}
					}
					if(min!=Double.MAX_VALUE){
						z=product(k_min,z);
					}
					else{
						break outer;
					}
					
					
				}
			}
			else{
				if(mod==1||mod==2){
					z=new double[]{0.01,0.01};
				}
				else if(mod==3||mod==4){
					z= new double[]{0.01,0.01,0.01,0.01};
				}
				else if(mod==5||mod==6){
					z= new double[]{0.01,0.01,0.01,0.01,0.01};
				}
			}
			
		
			x.add(diff(x.get(x.size()-1),z));//add scaling factor
			g.add(gradient(x.get(x.size()-1), mod, histo));
//			System.out.println((iter+1)+"	"+ou(x.get(x.size()-1))+"	"+model(histo,x.get(x.size()-1),mod));
			
			double val=model(histo,x.get(x.size()-1),mod);
			if(val<value_best_solution&&val>0&&in_standard_bounds(x.get(x.size()-1),standard_bounds[mod-1])){
				value_best_solution=val;
				best_solution=clone(x.get(x.size()-1));
			}
			
		}
		
		
		
		return concat(best_solution,value_best_solution);//concat(x.get(x.size()-1),model(histo,x.get(x.size()-1),mod));
	}
	
	
	/*
	public static double [] minimize_neg_ln_L(double[] start, int mod, double[][] bounds, ArrayList<int[]> histo ,int iterations, double step){//
		//bounds are currently ignored
		
		//double[] x=start;
		//double[] g=gradient(x);
		ArrayList<double[]> g=new ArrayList<double[]>();//stores the last m+1 entries
		ArrayList<double[]> x=new ArrayList<double[]>();//stores the last m+1 entries
		
		x.add(start);
		g.add(gradient(x.get(x.size()-1),mod, histo));
//		System.out.println((0)+"	"+ou(x.get(x.size()-1))+"	"+model(histo,x.get(x.size()-1),mod));
		
		double value_best_solution=model(histo,start,mod);
		double[] best_solution=clone(start);
		
		outer:
		for (int iter=0;iter<iterations;iter++){
			double[] z=new double[0];
			//double scaling=1;
			if(x.size()>1){
				double[] q=clone(g.get(g.size()-1));
				double [] alpha=new double[x.size()-1];
				for (int i=1;i<x.size();i++){
					double[] y=diff(g.get(x.size()-i),g.get(x.size()-i-1));
					double[] s=diff(x.get(x.size()-i),x.get(x.size()-i-1));
					double rho=1/product_tv_v(y,s);
					alpha[x.size()-i-1]=rho*product_tv_v(s,q);
					q=diff(q,product(alpha[i-1],y));
				}
				
				
				double[] y=diff(g.get(g.size()-1),g.get(g.size()-2));
				double[] s=diff(x.get(x.size()-1),x.get(x.size()-2));
				
				//double scaling=product_tv_v(s,y)/sq_norm(y);
				double[][] h_k_0=product(1.0/sq_norm(y),product_v_tv(y,s));
				z=product(h_k_0,q);
				for (int i=1;i<x.size();i++){
					y=diff(g.get(i),g.get(i-1));
					s=diff(x.get(i),x.get(i-1));
					double rho=1/product_tv_v(y,s);
					double beta=rho*product_tv_v(y,z);
					z=add(z,product(alpha[i-1]-beta,s));
				}
				//System.out.println("z: "+ou(z));
				if(is_nan(z)){
					break outer;
				}
				
				z=product(1/Math.sqrt(sq_norm(z)),z);
				
				double min=Double.MAX_VALUE;
				double k_min=model(histo,x.get(x.size()-1),mod);
				for (double k=0.01*step;k<=1*step;k+=0.01*step){
					double m=model(histo,diff(x.get(x.size()-1),product(k,z)),mod);
					if(m<min&&in_bound(diff(x.get(x.size()-1),product(k,z)),bounds)){
						min=m;
						k_min=k;
					}
				}
				if(min!=Double.MAX_VALUE){
					z=product(k_min,z);
				}
				else{
					min=Double.MAX_VALUE;
					k_min=model(histo,x.get(x.size()-1),mod);
					for (double k=0.0001*step;k<=0.01*step;k+=0.0001*step){
						double m=model(histo,diff(x.get(x.size()-1),product(k,z)),mod);
						if(m<min&&in_bound(diff(x.get(x.size()-1),product(k,z)),bounds)){
							min=m;
							k_min=k;
						}
					}
					if(min!=Double.MAX_VALUE){
						z=product(k_min,z);
					}
					else{
						break outer;
					}
					
					
				}
			}
			else{
				if(mod==1||mod==2){
					z=new double[]{0.01,0.01};
				}
				else if(mod==3||mod==4){
					z= new double[]{0.01,0.01,0.01,0.01};
				}
				else if(mod==5||mod==6){
					z= new double[]{0.01,0.01,0.01,0.01,0.01};
				}
			}
			
		
			x.add(diff(x.get(x.size()-1),z));//add scaling factor
			g.add(gradient(x.get(x.size()-1), mod, histo));
//			System.out.println((iter+1)+"	"+ou(x.get(x.size()-1))+"	"+model(histo,x.get(x.size()-1),mod));
			
			double val=model(histo,x.get(x.size()-1),mod);
			if(val<value_best_solution){
				value_best_solution=val;
				best_solution=clone(x.get(x.size()-1));
			}
			
		}
		
		
		
		return concat(best_solution,value_best_solution);//concat(x.get(x.size()-1),model(histo,x.get(x.size()-1),mod));
	}
	*/
	
	public static boolean is_nan(double[] a){
		for (int i=0;i<a.length;i++){
			if(Double.isNaN(a[i])){
				return true;
			}
		}
		return false;
	}
	
	
	public static boolean in_standard_bounds(double [] x, double[][] bounds){
		for (int i=0;i<bounds.length;i++){
			if(x[i]<=bounds[i][0]||bounds[i][1]<x[i]){
				return false;
			}
		}
		return true;
	}
	
	
	
	public static boolean in_bound(double [] x, double[][] bounds){
		for (int i=0;i<x.length;i++){
			if(x[i]<bounds[i][0]||bounds[i][1]<x[i]){
				return false;
			}
		}
		return true;
	}
	
	public static String ou(double[] a){
		String s=""+a[0];
		for (int i=1;i<a.length;i++){
			s=s+", "+a[i];
		}
		return s;
	}
	
	public static double[] concat(double [] x, double y){
		double [] xx=new double[x.length+1];
		for (int i=0;i<x.length;i++){
			xx[i]=x[i];
		}
		xx[xx.length-1]=y;
		return xx;
	}
	
	
	public static double[] gradient(double[] v, int mode, ArrayList<int[]> histo){
		double delta=0.0001;
		if(mode==1){
			return gradient1(v,histo,delta);
		}
		else if(mode==2){
			return gradient2(v,histo,delta);
		} 
		else if(mode==3){
			return gradient3(v,histo,delta);
		} 
		else if(mode==4){
			return gradient4(v,histo,0.0000001);
		} 
		else if(mode==5){
			return gradient5(v,histo,0.0000001);
		}
		else if(mode==6){
			return gradient6(v,histo,delta);
		} 
		else{
			return new double[0];
		}
	}
	
	public static double[] gradient1(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model1(histo,v[0]+delta,v[1])-model1(histo,v[0]-delta,v[1]))/(2*delta),(model1(histo,v[0],v[1]+delta)-model1(histo,v[0],v[1]-delta))/(2*delta)};
	}
	
	public static double[] gradient2(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model2(histo,v[0]+delta,v[1])-model2(histo,v[0]-delta,v[1]))/(2*delta),(model2(histo,v[0],v[1]+delta)-model2(histo,v[0],v[1]-delta))/(2*delta)};
	}
	
	public static double[] gradient3(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model3(histo,v[0]+delta,v[1],v[2],v[3])-model3(histo,v[0]-delta,v[1],v[2],v[3]))/(2*delta),(model3(histo,v[0],v[1]+delta,v[2],v[3])-model3(histo,v[0],v[1]-delta,v[2],v[3]))/(2*delta),(model3(histo,v[0],v[1],v[2]+delta,v[3])-model3(histo,v[0],v[1],v[2]-delta,v[3]))/(2*delta),(model3(histo,v[0],v[1],v[2],v[3]+delta)-model3(histo,v[0],v[1],v[2],v[3]-delta))/(2*delta)};
	}
	
	public static double[] gradient4(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model4(histo,v[0]+delta,v[1],v[2],v[3])-model4(histo,v[0]-delta,v[1],v[2],v[3]))/(2*delta),(model4(histo,v[0],v[1]+delta,v[2],v[3])-model4(histo,v[0],v[1]-delta,v[2],v[3]))/(2*delta),(model4(histo,v[0],v[1],v[2]+delta,v[3])-model4(histo,v[0],v[1],v[2]-delta,v[3]))/(2*delta),(model4(histo,v[0],v[1],v[2],v[3]+delta)-model4(histo,v[0],v[1],v[2],v[3]-delta))/(2*delta)};
	}
	
	public static double[] gradient5(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model5(histo,v[0]+delta,v[1],v[2],v[3],v[4])-model5(histo,v[0]-delta,v[1],v[2],v[3],v[4]))/(2*delta),(model5(histo,v[0],v[1]+delta,v[2],v[3],v[4])-model5(histo,v[0],v[1]-delta,v[2],v[3],v[4]))/(2*delta),(model5(histo,v[0],v[1],v[2]+delta,v[3],v[4])-model5(histo,v[0],v[1],v[2]-delta,v[3],v[4]))/(2*delta),(model5(histo,v[0],v[1],v[2],v[3]+delta,v[4])-model5(histo,v[0],v[1],v[2],v[3]-delta,v[4]))/(2*delta),(model5(histo,v[0],v[1],v[2],v[3],v[4]+delta)-model5(histo,v[0],v[1],v[2],v[3],v[4]-delta))/(2*delta)};
	}
	
	public static double[] gradient6(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model6(histo,v[0]+delta,v[1],v[2],v[3],v[4])-model6(histo,v[0]-delta,v[1],v[2],v[3],v[4]))/(2*delta),(model6(histo,v[0],v[1]+delta,v[2],v[3],v[4])-model6(histo,v[0],v[1]-delta,v[2],v[3],v[4]))/(2*delta),(model6(histo,v[0],v[1],v[2]+delta,v[3],v[4])-model6(histo,v[0],v[1],v[2]-delta,v[3],v[4]))/(2*delta),(model6(histo,v[0],v[1],v[2],v[3]+delta,v[4])-model6(histo,v[0],v[1],v[2],v[3]-delta,v[4]))/(2*delta),(model6(histo,v[0],v[1],v[2],v[3],v[4]+delta)-model6(histo,v[0],v[1],v[2],v[3],v[4]-delta))/(2*delta)};
	}
	
	
	public static double sq_norm(double[] v){
		double sum=0;
		for (int i=0;i<v.length;i++){
			sum+=v[i]*v[i];
		}
		return sum;
	}
	
	public static double[] clone(double[] a){
		double[] b=new double[a.length];
		for (int i=0;i<a.length;i++){
			b[i]=a[i];
		}
		return b;
	}
	
	public static double[] diff(double[] v, double[] w){
		double[] diff=new double[v.length];
		for (int i=0;i<v.length;i++){
			diff[i]=v[i]-w[i];
		}
		return diff;
	}
	
	public static double[] add(double[] v, double[] w){
		double[] add=new double[v.length];
		for (int i=0;i<v.length;i++){
			add[i]=v[i]+w[i];
		}
		return add;
	}
	
	public static double[] product(double[][] m, double[] v){
		double[] w=new double[m.length];
		for (int i=0;i<m.length;i++){
			for (int j=0;j<m[i].length;j++){
				w[i]+=m[i][j]*v[j];
			}
		}
		return w;
	}
	
	public static double[][] product_v_tv(double[] v, double[] w){
		double [][] m=new double[v.length][w.length];
		for (int i=0;i<v.length;i++){
			for (int j=0;j<w.length;j++){
				m[i][j]=v[i]*w[j];
			}
		}
		return m;
	}
	
	public static double product_tv_v(double[] v, double[] w){
		double m=0;
		for (int i=0;i<v.length;i++){
			for (int j=0;j<w.length;j++){
				m+=v[i]*w[j];
			}
		}
		return m;
	}

	public static double[][] product(double a, double[][] m){
		double[][] n=new double[m.length][m[0].length];
		for (int i=0;i<m.length;i++){
			for (int j=0;j<m[i].length;j++){
				n[i][j]=a*m[i][j];
			}
		}
		return n;
	}
	
	public static double[] product(double a, double[] v){
		double[] w=new double[v.length];
		for (int i=0;i<v.length;i++){
			w[i]=a*v[i];
		}
		return w;
	}
	
	
	public static ArrayList<int[]> histogram(ArrayList<Integer> counts){
		ArrayList<int[]> histogram=new ArrayList<int[]>();
		for(int i=0;i<counts.size();i++){
			int ii=index(counts.get(i).intValue(),histogram);
			if(ii==-1){
				histogram.add(new int[]{counts.get(i),1});
			}
			else{
				histogram.get(ii)[1]++;
			}
		}
		return histogram;
	}
	
	public static int index(int a, ArrayList<int[]> b){
		for (int i=0;i<b.size();i++){
			if(b.get(i)[0]==a){
				return i;
			}
		}
		return -1;
	}
	
	public static boolean is_olfactory(String s){
		
		return s.length()>=3&&s.substring(0,2).equals("OR")&&isInteger(s.substring(2,3));
	}
	
	public static boolean isInteger(String s){
		try{
			int a=Integer.parseInt(s);
			return true;
		}
		catch(Exception e){
			return false;
		}
	}
	
	
	//public static double[] minimize_neg_ln_L(double[] start, int mod_C, double[][] bounds){
	//	return new double[0];
	//}
	
	public static double[][] minimize_parallel(int rep_no, ArrayList<int[]> histo, int no_cpu, int iterations){
		ArrayList<double[]>[] solutions=new ArrayList[6];
		for (int mod=1;mod<=6;mod++){
			solutions[mod-1]=minimize_parallel(mod,rep_no,histo,no_cpu,iterations);
		}
		
		double[][][] standard_bounds={{{0,100},{0,100}},{{0,100},{0,100}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}}};
		
		SubthreadImprove[][] subthreads=new SubthreadImprove[solutions.length][];
		for (int i=0;i<solutions.length;i++){
			subthreads[i]=new SubthreadImprove[solutions[i].size()];
			for (int j=0;j<solutions[i].size();j++){
				subthreads[i][j]=new SubthreadImprove(i+1,  histo, solutions[i].get(j), standard_bounds[i],  iterations);
			}
		}
		execute(subthreads,no_cpu);
		
		double[][] solutions2=new double[6][];
		for (int i=0;i<subthreads.length;i++){
			double min=Double.MAX_VALUE;
			int j_min=-1;
			for (int j=0;j<subthreads[i].length;j++){
				if(subthreads[i][j].solution[subthreads[i][j].solution.length-1]<min){
					min=subthreads[i][j].solution[subthreads[i][j].solution.length-1];
					j_min=j;
				}
			}
			if(j_min!=-1){
				solutions2[i]=subthreads[i][j_min].solution;	
			}
			else{
				solutions2[i]=new double[0];
			}
			
		}
		return solutions2;
	}
	
	
	public static ArrayList<double[]> minimize_parallel(int mod_C, int rep_no, ArrayList<int[]> histo, int no_cpu, int iterations){
		double[][][] standard_bounds={{{0,100},{0,100}},{{0,100},{0,100}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}}};
		
		ArrayList<double[]> solutions=new ArrayList<double[]>();
		if(mod_C==1){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2)};
			
			//cur_min_res=new double[]{0,0,Double.MAX_VALUE};
			
			Subthread[] threads=new Subthread[rep_no];
			for (int i=0;i<rep_no;i++){
				threads[i] = new Subthread(new double[]{uniform(0.02,10), uniform(0.02,10)}, 1, new double[][]{{low_b[0],up_b[0]}, {low_b[1],up_b[1]}}, histo,iterations);	
			}
			execute(threads,no_cpu);
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i].solution,standard_bounds[0])){
					//System.out.println((i+1)+"	"+ou(p_res));
					if (threads[i].solution[2]>0 && threads[i].solution[2]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i].solution[0],threads[i].solution[1],threads[i].solution[2]});	
					}
				}
			}
			
			
		}
		else if(mod_C==2){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2)};
			//cur_min_res=new double[]{0,0,Double.MAX_VALUE};
			
			Subthread[] threads=new Subthread[rep_no];
			for (int i=0;i<rep_no;i++){
				threads[i]=new Subthread(new double[]{uniform(0.02,10), uniform(0.02,10)},2, new double[][]{{low_b[0],up_b[0]}, {low_b[1],up_b[1]}}, histo,iterations);
			}
			execute(threads,no_cpu);
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i].solution,standard_bounds[1])){
					//System.out.println((i+1)+"	"+ou(p_res));
					if (threads[i].solution[2]>0 && threads[i].solution[2]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i].solution[0],threads[i].solution[1],threads[i].solution[2]});	
					}
				}
			}
			
		}
		else if(mod_C==3){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			//cur_min_res=new double[]{0,0,0,0,Double.MAX_VALUE};
			
			Subthread[] threads=new Subthread[rep_no];
			for (int i=0;i<rep_no;i++){
				threads[i]=new Subthread(new double[]{uniform(0.02,10), uniform(0.02,10), uniform(0.02,10), uniform(2*Math.pow(10, -5),0.95)},  3, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],0.9999}},histo,iterations);
			}
			execute(threads,no_cpu);
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i].solution,standard_bounds[2])){
					if (threads[i].solution[4]>0 && threads[i].solution[4]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i].solution[0],threads[i].solution[1],threads[i].solution[2],threads[i].solution[3],threads[i].solution[4]});
					}	
				}	
			}
			
		}
		else if(mod_C==4){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			//cur_min_res=new double[]{0,0,0,0,Double.MAX_VALUE};
			
			Subthread[] threads=new Subthread[rep_no];
			for (int i=0;i<rep_no;i++){
				threads[i]=new Subthread(new double[]{uniform(0.02,10), uniform(0.02,10), uniform(0.02,10), uniform(2*Math.pow(10,-5),0.95)}, 4, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],0.9999}},histo,iterations);
				
			}
			execute(threads,no_cpu);
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i].solution,standard_bounds[3])){
					if (threads[i].solution[4]>0 && threads[i].solution[4]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i].solution[0],threads[i].solution[1],threads[i].solution[2],threads[i].solution[3],threads[i].solution[4]});
					}	
				}
			}
			
		}
		else if(mod_C==5){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			//cur_min_res=new double[]{0,0,0,0,0,Double.MAX_VALUE};
			
			Subthread[] threads=new Subthread[rep_no];
			for (int i=0;i<rep_no;i++){
				threads[i]=new Subthread(new double[]{uniform(0.02,10), uniform(0.02,5), uniform(0.02,10), uniform(0.02,10.),uniform(2*Math.pow(10, -5),0.95)},  5, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],up_b[3]},{low_b[4],0.9999}},histo,iterations);
			}
			execute(threads,no_cpu);
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i].solution,standard_bounds[4])){
					if (threads[i].solution[5]>0 && threads[i].solution[5]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i].solution[0],threads[i].solution[1],threads[i].solution[2],threads[i].solution[3],threads[i].solution[4],threads[i].solution[5]});
					}
				}
			}
			
		}	
		else if(mod_C==6){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			//cur_min_res=new double[]{0,0,0,0,0,Double.MAX_VALUE};
			
			Subthread[] threads=new Subthread[rep_no];
			for (int i=0;i<rep_no;i++){
				threads[i]=new Subthread(new double[]{uniform(0.02,10), uniform(0.02,5), uniform(0.02,10), uniform(0.02,10.),uniform(2*Math.pow(10, -5),0.95)},  6, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],up_b[3]},{low_b[4],0.9999}},histo,iterations);
			}
			execute(threads, no_cpu);
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i].solution,standard_bounds[5])){
					if (threads[i].solution[5]>0 && threads[i].solution[5]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i].solution[0],threads[i].solution[1],threads[i].solution[2],threads[i].solution[3],threads[i].solution[4],threads[i].solution[5]});
					}
				}
			}
		}
		Comparator<double[]> comp=(double[] a, double[] b)->{
			return new Double(a[a.length-1]).compareTo(b[b.length-1]);
		};
		Collections.sort(solutions,comp);
		
		ArrayList<double[]> solutions_final=new ArrayList<double[]>();
		for (int i=0;i<Math.min(5, solutions.size());i++){
			solutions_final.add(solutions.get(i));
		}	
		return solutions_final;
	}
	
	public static void execute(Subthread[] threads, int no_cpu){

		int no_undone=threads.length;
		int no_running=0;
		while(no_undone>0){
			no_undone=0;
			no_running=0;
			for (int i=0;i<threads.length;i++){
				if(!threads[i].done){
					no_undone++;
				}
				if(threads[i].running){
					no_running++;
				}
			}
			if(no_undone==0){
				break;
			}
			if(no_running<no_cpu){
				for (int i=0;i<threads.length;i++){
					if(!threads[i].done&&!threads[i].running){
						threads[i].start();
						no_running++;
					}
					if(no_running>=no_cpu){
						break;
					}
				}
			}
			try{
				Thread.sleep(500);
			}
			catch(Exception e){
				System.out.println(e);
			}
		}
	}
	
	
	public static void execute(SubthreadImprove[][] threads, int no_cpu){

		int no_undone=threads.length;
		int no_running=0;
		while(no_undone>0){
			no_undone=0;
			no_running=0;
			for (int i=0;i<threads.length;i++){
				for (int j=0;j<threads[i].length;j++){
					if(!threads[i][j].done){
						no_undone++;
					}
					if(threads[i][j].running){
						no_running++;
					}
				}
				
			}
			if(no_undone==0){
				break;
			}
			if(no_running<no_cpu){
				for (int i=0;i<threads.length;i++){
					for (int j=0;j<threads[i].length;j++){
						if(!threads[i][j].done&&!threads[i][j].running){
							threads[i][j].start();
							no_running++;
						}
						if(no_running>=no_cpu){
							break;
						}
					}
				}
			}
			try{
				Thread.sleep(500);
			}
			catch(Exception e){
				System.out.println(e);
			}
		}
	}
	
	private static class SubthreadImprove extends Thread{
		volatile double[] solution=new double[0];
		volatile boolean done=false;
		volatile boolean running=false;
		int mod_C=-1;
		ArrayList<int[]> histogram=new ArrayList<int[]>();
		double[] start=new double[0];
		double[][] bounds=new double[0][0];
		int iterations=0;
		
		public SubthreadImprove(int mod_C,  ArrayList<int[]> histogram, double[] start, double[][] bounds, int iterations){
			this.mod_C=mod_C;
			this.histogram=histogram;
			this.start=start;
			this.bounds=bounds;
			this.iterations=iterations;
		}
		
		public void run(){
			running=true;
			done=false;
			solution=improve(mod_C, histogram, start, bounds, iterations);
			running=false;
			done=true;
		}
	}
	
	private static class Subthread extends Thread{
		volatile double[] solution=new double[0];
		volatile boolean done=false;
		volatile boolean running=false;
		double[] start=new double[0];
		int mod=-1;
		double[][] bounds=new double[0][0];
		ArrayList<int[]> histo=new ArrayList<int[]>();
		int iterations=0;
		public Subthread(double[] start, int mod, double[][] bounds, ArrayList<int[]> histo, int iterations){
			this.start=start;
			this.mod=mod;
			this.bounds=bounds;
			this.histo=histo;
			this.iterations=iterations;
		}
		public void run(){
			running=true;
			done=false;
			solution=minimize_neg_ln_L( start,  mod, bounds,  histo, iterations);
			running=false;
			done=true;
		}
	}
	
	
	
	public static double[] minimize(int mod_C, int rep_no, ArrayList<int[]> histo, int iterations){
		double[] cur_min_res=new double[0];
		if(mod_C==1){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2)};
			double[][] standard_bounds={{0,100},{0,100}};
			cur_min_res=new double[]{0,0,Double.MAX_VALUE};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,10)}, 1, new double[][]{{low_b[0],up_b[0]}, {low_b[1],up_b[1]}}, histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[2]>0 && p_res[2]<cur_min_res[2]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2]};	
					}
				}		
			}
			
		}
		else if(mod_C==2){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2)};
			cur_min_res=new double[]{0,0,Double.MAX_VALUE};
			double[][] standard_bounds={{0,100},{0,100}};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,10)},2, new double[][]{{low_b[0],up_b[0]}, {low_b[1],up_b[1]}}, histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[2]>0 && p_res[2]<cur_min_res[2]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2]};	
					}
				}
				
			}
			
		}
		else if(mod_C==3){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			cur_min_res=new double[]{0,0,0,0,Double.MAX_VALUE};
			double[][] standard_bounds={{0,100},{0,100},{0,100},{0,1}};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,10), uniform(0.02,10), uniform(2*Math.pow(10, -5),0.95)},  3, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],0.9999}},histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[4]>0 && p_res[4]<cur_min_res[4]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2],p_res[3],p_res[4]};
					}
				}
				
			}
		}
		else if(mod_C==4){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			cur_min_res=new double[]{0,0,0,0,Double.MAX_VALUE};
			double[][] standard_bounds={{0,100},{0,100},{0,100},{0,1}};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,10), uniform(0.02,10), uniform(2*Math.pow(10,-5),0.95)}, 4, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],0.9999}},histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[4]>0 && p_res[4]<cur_min_res[4]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2],p_res[3],p_res[4]};
					}
				}
				
			}
		}
		else if(mod_C==5){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			cur_min_res=new double[]{0,0,0,0,0,Double.MAX_VALUE};
			double[][] standard_bounds={{0,100},{0,100},{0,100},{0,100},{0,1}};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,5), uniform(0.02,10), uniform(0.02,10.),uniform(2*Math.pow(10, -5),0.95)},  5, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],up_b[3]},{low_b[4],0.9999}},histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[5]>0 && p_res[5]<cur_min_res[5]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2],p_res[3],p_res[4],p_res[5]};
					}
				}
				
			}
		}	
		else if(mod_C==6){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			cur_min_res=new double[]{0,0,0,0,0,Double.MAX_VALUE};
			double[][] standard_bounds={{0,100},{0,100},{0,100},{0,100},{0,1}};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,5), uniform(0.02,10), uniform(0.02,10.),uniform(2*Math.pow(10, -5),0.95)},  6, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],up_b[3]},{low_b[4],0.9999}},histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[5]>0 && p_res[5]<cur_min_res[5]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2],p_res[3],p_res[4],p_res[5]};
						
					}
				}
				
			}
		}
		return cur_min_res;
	}
	
	
	
	public static double model1(int [] s, double a, double b){
		double sum=0;
		for (int i=0;i<s.length;i++){
			sum += (s[i]*Math.log(b) + (-s[i]-a)*Math.log(1 + b) + Gamma.logGamma(s[i] + a) -  Gamma.logGamma(s[i]+1) -  Gamma.logGamma(a));
		}
		return -sum;
	}
	
	public static double model(ArrayList<int[]> histo, double[] param, int mod){
		if(mod==1){
			return model1(histo,param[0],param[1]);
		}
		else if(mod==2){
			return model2(histo,param[0],param[1]);
		}
		else if(mod==3){
			return model3(histo,param[0],param[1],param[2],param[3]);
		}
		else if(mod==4){
			return model4(histo,param[0],param[1],param[2],param[3]);
		}
		else if(mod==5){
			return model5(histo,param[0],param[1],param[2],param[3],param[4]);
		}
		else if(mod==6){
			return model6(histo,param[0],param[1],param[2],param[3],param[4]);
		}
		else{
			return 0;
		}
	}
	
	public static double model1(ArrayList<int[]> h, double a, double b){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x= (h.get(i)[0]*Math.log(b) + (-h.get(i)[0]-a)*Math.log(1 + b) + Gamma.logGamma(h.get(i)[0] + a) -  Gamma.logGamma(h.get(i)[0]+1) -  Gamma.logGamma(a));
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
		}
		return -sum;
	}
	
	public static double model2(int [] s, double a, double b){
		double sum=0;
		for (int i=0;i<s.length;i++){
			double proxy=besselk(-s[i] + a, 2*Math.sqrt(b));
			sum+= Math.log(2) + ((s[i] + a)/2.0)*Math.log(b) + Math.log(proxy) - Gamma.logGamma(s[i]+1) - Gamma.logGamma(a);
		}
		return -sum;		
	}
	public static double model2(ArrayList<int[]> h, double a, double b){//TODO: bessel
		double sum=0;
		for (int i=0;i<h.size();i++){
			double proxy=besselk(-h.get(i)[0] + a, 2*Math.sqrt(b));
			//System.out.println(h.get(i)[0]+"	"+proxy);
			if(proxy>0&&Double.isFinite(proxy)){
				double x= Math.log(2) + ((h.get(i)[0] + a)/2.0)*Math.log(b) + Math.log(proxy) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(a);
				if(Double.isFinite(x)){
					sum+=x*h.get(i)[1];
				}
				
			}
		}
		return -sum;		
	}
	
	public static double model3(int [] s, double a, double b, double t, double w){
		double sum=0;
		for (int i=0;i<s.length;i++){
			sum += Math.log( Math.exp( Math.log(w * t) + (-1-s[i])*Math.log(1 + t) ) + Math.exp( Math.log(1-w) + s[i]*Math.log(b) + (-s[i]-a)*Math.log(1 + b) + Gamma.logGamma(s[i] + a) - Gamma.logGamma(s[i]+1) - Gamma.logGamma(a) ) );
		}
		return -sum;
	}
	public static double model3(ArrayList<int[]> h, double a, double b, double t, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x= Math.log( Math.exp( Math.log(w * t) + (-1-h.get(i)[0])*Math.log(1 + t) ) + Math.exp( Math.log(1-w) + h.get(i)[0]*Math.log(b) + (-h.get(i)[0]-a)*Math.log(1 + b) + Gamma.logGamma(h.get(i)[0] + a) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(a) ) );
			//System.out.println(h.get(i)[0]+"	"+x);
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
		}
		return -sum;
	}
	
	public static double model4(int [] s, double a, double b, double t, double w){
		double sum=0;
		for (int i=0;i<s.length;i++){
			if (s[i]>25){
				sum += Math.log( Math.exp( Math.log(w * t) + (-1 - s[i])*Math.log(1 + t) ) + Math.exp( Math.log(1-w) + Math.log(2) + ((s[i] + a)/2.0)*Math.log(b) + Math.log(besselk(-s[i] + a, 2*Math.sqrt(b))) - Gamma.logGamma(s[i] + 1) - Gamma.logGamma(a) ) );	
			}
			else{
				sum += Math.log( Math.exp( Math.log(w * t) + (-1 - s[i])*Math.log(1 + t) ) + Math.exp( Math.log(1-w) + Math.log(2) + ((s[i] + a)/2.0)*Math.log(b) + Math.log(kv(-s[i] + a, 2*Math.sqrt(b))) - Gamma.logGamma(s[i] + 1) - Gamma.logGamma(a) ) );
						
			}
		}
		return -sum;
	}
	
	public static double model4(ArrayList<int[]> h, double a, double b, double t, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			if (h.get(i)[0]>25){
				double x= Math.log( Math.exp( Math.log(w * t) + (-1 - h.get(i)[0])*Math.log(1 + t) ) + Math.exp( Math.log(1-w) + Math.log(2) + ((h.get(i)[0] + a)/2.0)*Math.log(b) + Math.log(besselk(-h.get(i)[0] + a, 2*Math.sqrt(b))) - Gamma.logGamma(h.get(i)[0] + 1) - Gamma.logGamma(a) ) );	
				if(Double.isFinite(x)){
					sum+=x*h.get(i)[1];
				}
			}
			else{
				double x= Math.log( Math.exp( Math.log(w * t) + (-1 - h.get(i)[0])*Math.log(1 + t) ) + Math.exp( Math.log(1-w) + Math.log(2) + ((h.get(i)[0] + a)/2.0)*Math.log(b) + Math.log(kv(-h.get(i)[0] + a, 2*Math.sqrt(b))) - Gamma.logGamma(h.get(i)[0] + 1) - Gamma.logGamma(a) ) );
				if(Double.isFinite(x)){
					sum+=x*h.get(i)[1];		
				}
			}
		}
		return -sum;
	}
	
	public static double model5(int [] s, double a, double b, double g, double d, double w){
		double sum=0;
		for (int i=0;i<s.length;i++){
			sum += Math.log( Math.exp( Math.log(w) + s[i]*Math.log(b) + (-s[i]-a)*Math.log(1 + b) + Gamma.logGamma(s[i] + a) - Gamma.logGamma(s[i]+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1-w) + s[i]*Math.log(d) + (-s[i]-g)*Math.log(1 + d) + Gamma.logGamma(s[i] + g) - Gamma.logGamma(s[i]+1) - Gamma.logGamma(g) ) );
		}
		return -sum;
	}
	
	public static double model5(ArrayList<int[]> h, double a, double b, double g, double d, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x= Math.log( Math.exp( Math.log(w) + h.get(i)[0]*Math.log(b) + (-h.get(i)[0]-a)*Math.log(1 + b) + Gamma.logGamma(h.get(i)[0] + a) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1-w) + h.get(i)[0]*Math.log(d) + (-h.get(i)[0]-g)*Math.log(1 + d) + Gamma.logGamma(h.get(i)[0] + g) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(g) ) );
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
		}
		return -sum;
	}
	
	public static double model6(int [] s, double a, double b, double g, double d, double w){
		double sum=0;
		for (int i=0;i<s.length;i++){
			if (s[i]>25){
				sum += Math.log( (w * Math.exp(s[i]*Math.log(b) + (-s[i]-a)*Math.log(1 + b) + Gamma.logGamma(s[i] + a) - Gamma.logGamma(s[i]+1) - Gamma.logGamma(a))) + ((1-w) * Math.exp( Math.log(2) + ((s[i] + g)/2.0)*Math.log(d) + Math.log(besselk(-s[i] + g, 2*Math.sqrt(d))) - Gamma.logGamma(s[i] + 1) - Gamma.logGamma(g) ) ));		
			}
			else{
				sum += Math.log( (w * Math.pow(b,s[i]) * Math.pow(1 + b,-s[i]-a) * Math.exp(Gamma.logGamma(s[i] + a) - Gamma.logGamma(s[i]+1) - Gamma.logGamma(a))) + ((1-w) * Math.exp( Math.log(2) + ((s[i] + g)/2.0)*Math.log(d) + Math.log(kv(-s[i] + g, 2*Math.sqrt(d))) - Gamma.logGamma(s[i] + 1) - Gamma.logGamma(g) ) ));

			}
		}
		return -sum;
	}
	
	public static double model6(ArrayList<int[]> h, double a, double b, double g, double d, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			if (h.get(i)[0]>25){
				double x= Math.log( (w * Math.exp(h.get(i)[0]*Math.log(b) + (-h.get(i)[0]-a)*Math.log(1 + b) + Gamma.logGamma(h.get(i)[0] + a) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(a))) + ((1-w) * Math.exp( Math.log(2) + ((h.get(i)[0] + g)/2.0)*Math.log(d) + Math.log(besselk(-h.get(i)[0] + g, 2*Math.sqrt(d))) - Gamma.logGamma(h.get(i)[0] + 1) - Gamma.logGamma(g) ) ));		
				if(Double.isFinite(x)){
					sum+=x*h.get(i)[1];
				}
			}
			else{
				double x= Math.log( (w * Math.pow(b,h.get(i)[0]) * Math.pow(1 + b,-h.get(i)[0]-a) * Math.exp(Gamma.logGamma(h.get(i)[0] + a) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(a))) + ((1-w) * Math.exp( Math.log(2) + ((h.get(i)[0] + g)/2.0)*Math.log(d) + Math.log(kv(-h.get(i)[0] + g, 2*Math.sqrt(d))) - Gamma.logGamma(h.get(i)[0] + 1) - Gamma.logGamma(g) ) ));
				if(Double.isFinite(x)){
					sum+=x*h.get(i)[1];
				}
			}
		}
		return -sum;
	}
	
	public static double uniform (double a, double b){
		return a+Math.random()*(b-a);
	}
	public static double besselk(double n, double x){
		return Bessel.k(x, n, false);
	}
	public static double kv(double n, double x){
		return Bessel.k(x, n, false);
	}
}
