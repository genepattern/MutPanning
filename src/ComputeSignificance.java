/************************************************************           
 * MutPanning - Step 11										*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2018 									*   
 *															*   
 * License:		Public Domain								*   
 *															*   
 * Summary: This script computes the mutational 			*
 *  significance of each gene. The major component of this	*
 *  test is based on the combined null hypothesis model,	*
 *  which integrates mutation counts and sequence context	*
 *  around mutations. Further this method tests for			*
 *  local clustering of mutations in mutation hotspots		*
 *  (enhancing the sensitivity for oncogenes)and for 		*
 *  accumulation of protein damaging mutations (enhancing	*
 *  the sensitivity for tumor supressor genes). With the	*
 *  latter test insersions and deletions are considered in 	*
 *  our test statistics.									*
 *************************************************************/

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.special.Gamma;



public class ComputeSignificance {
	static String[] chr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};	
	static ArrayList<Gene>[] genes=new ArrayList[chr.length];
	
	static String entity="";
	static String file_genes="";
	static String file_count="";
	static String file_cbase="";
	static String file_destructive="";
	static String file_gene_count="";
	static String file_count_genes="";
	static String file_amino="";
	static String file_out="";
	
	static String[] exclude={
			"PPIAL4A",//too long
			"PPIAL4B",
			"PPIAL4C",
			"FAM75A5",
			"FAM75A7",
			"ANKRD20A3",
			"FOXD4L2",
			"FOXD4L4",
			"FAM21B",
			"ASMT",//appears in file twice
			"ASMTL",
			"CSF2RA",
			"DHRSX",
			"IL3RA",
			"IL3R",
			"IL9R",
			"SPRY3",
			"ZBED1"};
	
	/*
	 * argument0 : root file
	 * argument1: entity
	* argument2: path to Hg19 folder
	 */
	
	public static void main(String[] args){
		try{
			entity=args[1];
			file_genes=args[2]+”Exons_Hg19.txt";
			file_count=args[0]+"EntityCounts/EntityCounts";
			file_cbase=args[0]+"CBASE/Distributions/";
			file_destructive=args[0]+"CountDestructive/";
			file_gene_count=args[0]+"CBASE/Counts/Count";
			file_count_genes=args[0]+"CBASE/CountsChrwise/Count";
			file_amino=args[2]+”ASAnnotationHg19/ASAnnotation_chr";
			file_out=args[0]+"SignificanceRaw/Significance";
			
			if(!new File(args[0]+"SignificanceRaw/").exists()){
				new File(args[0]+"SignificanceRaw/").mkdir();
			}
			
			for (int i=0;i<genes.length;i++){
				genes[i]=new ArrayList<Gene>();
			}
			
			//read genes with coordinates
			FileInputStream in=new FileInputStream(file_genes);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			String s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				if(index(t[0],exclude)!=-1){
					continue;
				}
				int ii=index_gene(t[0],genes[Integer.parseInt(t[1])-1]);
				if(ii==-1){
					genes[Integer.parseInt(t[1])-1].add(new Gene(t[0],Integer.parseInt(t[2]),Integer.parseInt(t[3])));
				}
				else{
					genes[Integer.parseInt(t[1])-1].get(ii).coord.add(new int[]{Integer.parseInt(t[2]),Integer.parseInt(t[3])});
				}
			}
			input.close();
			
			//define external intervals of gene (needed in the subthread to see when the sequence of a gene is fully read
			//and computation of the significance values can start
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					genes[i].get(j).start=min(genes[i].get(j).coord);
					genes[i].get(j).end=max(genes[i].get(j).coord);
				}
			}
			
			
			//these are the distribution of the total mutation count computed for each gene
			//based on the background distribution of the synonymous mutations across the genome
			in=new FileInputStream(file_cbase+entity+".txt");
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int i1=-1;
			int i2=-1;
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				if(t.length==1){
					i1=-1;
					i2=-1;
					boolean f=false;
					outer:
					for (int i=0;i<genes.length;i++){
						for (int j=0;j<genes[i].size();j++){
							if(genes[i].get(j).name.equals(t[0])){
								//found[i][j]=true;
								i1=i;
								i2=j;
								f=true;
								break outer;
							}
						}
					}
					if(!f){
						String transform=transform(t[0]);//avoid the artifacts introduced by Excel 
						f=false;
						outer:
						for (int i=0;i<genes.length;i++){
							for (int j=0;j<genes[i].size();j++){
								if(genes[i].get(j).name.equals(transform)){
									//found[i][j]=true;
									i1=i;
									i2=j;
									f=true;
									break outer;
								}
							}
						}
					}
					
				}
				else if(i1!=-1&&i2!=-1){
					if(Integer.parseInt(t[0])!=genes[i1].get(i2).prob_nonsyn.size()){
						System.out.println("Nonsyn File Incompatible");
						System.out.println(genes[i1].get(i2).name);
						System.exit(0);
					}
					genes[i1].get(i2).prob_nonsyn.add(Double.parseDouble(t[1]));
				}
			}
			input.close();
			
			//normalize the distributions
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					double sum=0;
					for (int k=0;k<genes[i].get(j).prob_nonsyn.size();k++){
						sum+=genes[i].get(j).prob_nonsyn.get(k);
					}
					if(sum>1){
						for (int k=0;k<genes[i].get(j).prob_nonsyn.size();k++){
							genes[i].get(j).prob_nonsyn.set(k,genes[i].get(j).prob_nonsyn.get(k)/sum);
						}
					}
				}
			}
			
			//read the observed syn and nonsyn counts for each gene
			for (int i=0;i<chr.length;i++){
				in=new FileInputStream(file_count_genes+entity+"_Chr"+chr[i]+".txt");
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					int ii=index_gene(t[0],genes[i]);
					if(ii!=-1){
						genes[i].get(ii).cov=Double.parseDouble(t[1]);
						genes[i].get(ii).cov_syn=Double.parseDouble(t[2]);
						genes[i].get(ii).count=Integer.parseInt(t[3]);
						genes[i].get(ii).count_syn=Integer.parseInt(t[4]);
					}
					
				}
				input.close();
			}
			
			
			Comparator<Gene> comp=(Gene a, Gene b)->{
				return new Integer(a.start).compareTo(new Integer(b.start));
			};
			
			for (int i=0;i<genes.length;i++){
				Collections.sort(genes[i],comp);
			}
			
			
			//compute the significance for genes per chr in paralle
			Subthread[] threads=new Subthread[chr.length];
			for (int i=0;i<threads.length;i++){
				threads[i]=new Subthread();
				threads[i].c=i;
				threads[i].start();
				
			}

			//wait until all computations are done
			boolean all_done=false;
			do{
				Thread.sleep(3000);
				all_done=true;
				for (int i=0;i<threads.length;i++){
					if(!threads[i].done){
						all_done=false;
						break;
					}
				}
			}while(!all_done);
			
			
			//as the last component we further test for accumulations of destructive
			//mutations in genes (including insertions/deletions)
			
			//read no. destructive mutations
			int count_all=0;
			int count_destruct=0;
			in=new FileInputStream(file_destructive+entity+".txt");
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				int[] index=index(t[0],genes);
				if(index[0]!=-1){
					genes[index[0]].get(index[1]).count_all=Integer.parseInt(t[1]);
					genes[index[0]].get(index[1]).count_destruct=Integer.parseInt(t[2]);
					count_all+=Integer.parseInt(t[1]);
					count_destruct+=Integer.parseInt(t[2]);
				}
			}
			input.close();
			
			//test for accumulation of destructive mutations (just a very straightforward binomial distribution does the job)	
			double frac_destruct=(double)(count_destruct)/(double)(count_all);
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					BinomialDistribution dist=new BinomialDistribution(genes[i].get(j).count_all,frac_destruct);
					if(genes[i].get(j).count_all<2||genes[i].get(j).count_destruct==0){
						genes[i].get(j).sign_destruct=1-dist.cumulativeProbability(genes[i].get(j).count_destruct-1);
					}
					else{
						genes[i].get(j).sign_destruct=1-dist.cumulativeProbability(genes[i].get(j).count_destruct-1);
					}
				}
			}
			
			//combination of the 3 p-values using brown
			ArrayList<Double> product_syn=new ArrayList<Double>();
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					double prod=genes[i].get(j).sign_vector_syn*genes[i].get(j).sign_hotspot_syn;
					if(prod!=0&&prod!=1){
						product_syn.add(prod);
					}
				}
			}
			
			
			double avg_syn=0;
			for (int i=0;i<product_syn.size();i++){
				avg_syn+=(-2)*Math.log(product_syn.get(i));
			}
			avg_syn/=(double)(product_syn.size());
			
			double var_syn=0;
			for (int i=0;i<product_syn.size();i++){
				var_syn+=((-2)*Math.log(product_syn.get(i))-avg_syn)*((-2)*Math.log(product_syn.get(i))-avg_syn);
			}
			var_syn/=(double)(product_syn.size());
			
			//determine the right no. dimensions and scaling factor for the Brown method for the syn mutations
			double c1_syn=var_syn/(2*avg_syn);
			double k1_syn=2*avg_syn*avg_syn/var_syn;
			
			
			//Combine p-values using brown for the syn mutations
			ChiSquaredDistribution dist1_syn=new ChiSquaredDistribution(k1_syn);
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					double prod=genes[i].get(j).sign_vector_syn*genes[i].get(j).sign_hotspot_syn;
					if(prod==1){
						genes[i].get(j).sign_complete_syn=1;
					}
					else if(prod==0){
						genes[i].get(j).sign_complete_syn=0;
					}
					else{
						genes[i].get(j).sign_complete_syn=1-dist1_syn.cumulativeProbability((-2)*Math.log(prod)/c1_syn);
					}
				}
			}
			
			
			
			ArrayList<Double> product=new ArrayList<Double>();
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					double prod=genes[i].get(j).sign_combined*Math.min(genes[i].get(j).sign_destruct,genes[i].get(j).sign_hotspot);
					if(prod!=0&&prod!=1){
						product.add(prod);
					}
				}
			}
			
			double avg=0;
			for (int i=0;i<product.size();i++){
				avg+=(-2)*Math.log(product.get(i));
			}
			avg/=(double)(product.size());
			
			double var=0;
			for (int i=0;i<product.size();i++){
				var+=((-2)*Math.log(product.get(i))-avg)*((-2)*Math.log(product.get(i))-avg);
			}
			var/=(double)(product.size());
			
			
			//determine the right no. dimensions and scaling factor for the Brown method for the nonsyn mutations
			double c1=var/(2*avg);
			double k1=2*avg*avg/var;
			//Combine p-values using brown for the nonsyn mutations
			
			ChiSquaredDistribution dist1=new ChiSquaredDistribution(k1);
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					double prod=genes[i].get(j).sign_combined*Math.min(genes[i].get(j).sign_destruct,genes[i].get(j).sign_hotspot);
					if(prod==1){
						genes[i].get(j).sign_complete=1;
					}
					else if(prod==0){
						genes[i].get(j).sign_complete=0;
					}
					else{
						genes[i].get(j).sign_complete=1-dist1.cumulativeProbability((-2)*Math.log(prod)/c1);
					}
				}
			}
			
			//for FDR correction put all genes in 1 list and sort them by sign (both for syn sign and non syn sign)
			ArrayList<Gene> genes_all=new ArrayList<Gene>();
			for (int i=0;i<genes.length;i++){
				genes_all.addAll(genes[i]);
			}
			
			
			Comparator<Gene> comp_gene_syn=(Gene g1, Gene g2)->{
				return new Double(g1.sign_complete_syn).compareTo(new Double(g2.sign_complete_syn));
			};
			
			Collections.sort(genes_all,comp_gene_syn);
			for (int i=0;i<genes_all.size();i++){
				genes_all.get(i).fdr_syn=Math.min(1, (double)(genes_all.size())/(double)(i+1)*genes_all.get(i).sign_complete_syn);
			}
			
			
			Comparator<Gene> comp_gene=(Gene g1, Gene g2)->{
				return new Double(g1.sign_complete).compareTo(new Double(g2.sign_complete));
			};
			
			
			int NN=0;
			for (int i=0;i<genes_all.size();i++){
				if(genes_all.get(i).count_syn+genes_all.get(i).count+genes_all.get(i).count_all+genes_all.get(i).count_destruct>0){
					NN++;
				}
			}
			
			//Compute fdr
			Collections.sort(genes_all,comp_gene);
			for (int i=0;i<genes_all.size();i++){
				if(i+1<=NN){
					genes_all.get(i).fdr=Math.min(1, (double)(NN)/(double)(i+1)*genes_all.get(i).sign_complete);
				}
				else{
					//genes_all.get(i).fdr=1;
					break;
				}
				
			}
			
			//output of FDR values and 
			FileWriter out=new FileWriter(file_out+entity+".txt");
			BufferedWriter output= new BufferedWriter(out);
			
			output.write("Name	TargetSize	TargetSizeSyn	Count	CountSyn	SignVectorSyn	SignCumSyn	SignSeqSyn	SignVector	SignCount	SignSeq	SignDm	SignCum	SignCombined	FDRSyn	FDR");
			output.newLine();
			
			for (int i=0;i<genes_all.size();i++){
				output.write(genes_all.get(i).name+"	"+genes_all.get(i).cov+"	"+genes_all.get(i).cov_syn+"	"+genes_all.get(i).count+"	"+genes_all.get(i).count_syn+"	"+genes_all.get(i).sign_vector_syn+"	"+genes_all.get(i).sign_hotspot_syn+"	"+genes_all.get(i).sign_complete_syn+"	"+genes_all.get(i).sign_vector+"	"+genes_all.get(i).sign_cbase+"	"+genes_all.get(i).sign_combined+"	"+genes_all.get(i).sign_destruct+"	"+genes_all.get(i).sign_hotspot+"	"+genes_all.get(i).sign_complete+"	"+genes_all.get(i).fdr_syn+"	"+genes_all.get(i).fdr);
				output.newLine();
			}
			output.close();
			
		}
		catch(Exception e){
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			System.out.println(e);
		}
		
	}
	
	

	
	public static int[] index (String name, ArrayList<Gene>[] genes){
		for (int i=0;i<genes.length;i++){
			for (int j=0;j<genes[i].size();j++){
				if(genes[i].get(j).name.equals(name)){
					return new int[]{i,j};
				}
			}
		}
		return new int[]{-1,-1};
	}
	
	public static String transform(String s){
		if(s.contains("-Mar")){
			return "MARCH"+s.split("-")[0];
		}
		else if(s.contains("-Sep")){
			if(s.split("-")[0].equals("15")){
				return "SEP"+s.split("-")[0];
			}
			return "SEPT"+s.split("-")[0];
		}
		else if(s.contains("-Dec")){
			return "DEC"+s.split("-")[0];
		}
		else{
			return s;
		}
		
	}
	
	
	//an object to store the genomic coordinates and statistical test results
	private static class Gene{
		String name="";
		ArrayList<int[]> coord=new ArrayList<int[]>();
		int start=-1;
		int end=-1;
		ArrayList<Double> prob_nonsyn=new ArrayList<Double>();
		
		double cov=0;
		double cov_syn=0;
		int count=0;
		int count_syn=0;
		
		//double sign_syn=1;
		
		int count_all=0;
		int count_destruct=0;
		
		
		double sign_vector=1;
		double sign_cbase=1;
		double sign_combined=1;
		double sign_destruct=1;
		double sign_hotspot=1;
		double sign_complete=1;
		double sign_vector_syn=1;
		double sign_hotspot_syn=1;
		double sign_complete_syn=1;
		
		double fdr=1;
		double fdr_syn=1;
		public Gene(String name, int start, int end){
			this.name=name;
			coord.add(new int[]{start,end});
		}
		public boolean contains (int pos){
			for (int i=0;i<coord.size();i++){
				if(coord.get(i)[0]<=pos&&pos<=coord.get(i)[1]){
					return true;
				}
			}
			return false;
		}
	}
	
	public static int[] index(ArrayList<Double> a){
		//double sum=0;
		ArrayList<double[]> index=new ArrayList<double[]>();
		for (int i=0;i<a.size();i++){
			index.add(new double[]{a.get(i),i});
		}
		Comparator<double[]> comp=(double[] x, double[] y)->{
			return -new Double(x[0]).compareTo(new Double(y[0]));
		};
		Collections.sort(index,comp);
		int[] ii=new int[index.size()];
		for (int i=0;i<index.size();i++){
			ii[i]=(int)(index.get(i)[1]);
		}
		return ii;
	}
	
	public static boolean contains(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
	public static int index(String s, String[] t){
		for (int i=0;i<t.length;i++){
			if(t[i].equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	public static int min(ArrayList<int[]> a){
		int min=1000000000;
		for (int i=0;i<a.size();i++){
			if(a.get(i)[0]<min){
				min=a.get(i)[0];
			}
		}
		return min;
	}
	
	public static int max(ArrayList<int[]> a){
		int max=-1000000000;
		for (int i=0;i<a.size();i++){
			if(a.get(i)[1]>max){
				max=a.get(i)[1];
			}
		}
		return max;
	}
	
	public static int index(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return i;
			}
		}
		return -1;
	}
	public static int index_gene(String s, ArrayList<Gene> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).name.equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	
	//the thread that is responsible to compute the significance
	//in brief reads for each gene the mutation rates of its positions
	//computes the probability to observe the distribution under a multinomial 
	//distribution and derives the p-value by using Monte-Carlo
	private static class Subthread extends Thread{
		volatile boolean done=false;
		int c=-1;
		ThreadLocalRandom random_conc=null;
		public double sign(ArrayList<Double> lambda2, ArrayList<Integer> count2 ,int err){
			double sign_quick=sign_quick(lambda2,count2, err);
			if(sign_quick<0.05){
				return sign_long(lambda2,count2, err);
			}
			else{
				return sign_quick;
			}
			//return sign_quick;
		}
		
		public double sign(ArrayList<Double> lambda2, ArrayList<Integer> count2 ,ArrayList<Double> probability){//int err, 
			double sign_quick=sign_quick(lambda2,count2, probability);//err,
			if(sign_quick<0.05){
				double d=sign_long(lambda2,count2,probability);
				//System.out.println("sign	"+sign_quick+"	"+d);
				return d;//err
			}
			else{
				return sign_quick;
			}
		
		}
		
		public double signREVISED(int count, ArrayList<Double> lambda2, ArrayList<Integer> count2 ,ArrayList<Double> probability, int err){//int err, 
			double sign_quick=sign_quickREVISED(count, lambda2,count2, probability,err);//err,
			if(sign_quick<0.05){
				double d=sign_longREVISED(count, lambda2,count2,probability,err);
				//System.out.println("sign	"+sign_quick+"	"+d);
				return d;//err
			}
			else{
				return sign_quick;
			}
			//return sign_quick;
		}
		
		
		
		public int[] random_collision (int N, double[] cum){
			int[] count=new int[cum.length];
			for (int i=0;i<N;i++){
				
				double r=random_conc.nextDouble();//Math.rrandom();
				for (int j=0;j<cum.length;j++){
					if(j<cum.length-1){
						if(cum[j]<=r&&r<cum[j+1]){
							count[j]++;
							break;
						}
					}
					else{
						if(cum[cum.length-1]<=r){
							count[cum.length-1]++;
							break;
						}
					}

				}
			}
		
			int[] random_meta=new int[10];
			for (int i=0;i<count.length;i++){
				if(count[i]>0&&count[i]<=random_meta.length){
					random_meta[count[i]-1]++;
				}
			}
			return random_meta;
		}
		
		public double sign_long(ArrayList<Double> lambda2, ArrayList<Integer> count2, int err){
			
			
			int NN=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
			}
			if(lambda2.size()==0||NN==0){
				return 1;
			}
			
			double T=0;
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
					//System.out.println(Math.log(lambda2.get(i)));
				}
				if(count2.get(i)>0){
					T-=Gamma.logGamma(count2.get(i)+1);
				}
			}
			double offset=NN*Math.log((double)(1)/(double)(lambda2.size()));
			T+=offset;
			//System.out.println(T);
			
			//System.out.println("T: "+Math.exp(T));
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			//int NN=500;
			
			Collections.sort(lambda2,comp);
			double sum_lambda=0;
			for (int i=0;i<lambda2.size();i++){
				sum_lambda+=lambda2.get(i);
			}
			
			double sum_pow=0;
			{
				int k=2;
				for (int i=0;i<lambda2.size();i++){
					sum_pow+=Math.pow(lambda2.get(i),k);
				}
			}
		
			

			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			{
				int k=2;
				double sum_local=0;
				for (int i=0;i<lambda2.size();i++){
					sum_local+=Math.pow(lambda2.get(i),k);
					if(sum_local/sum_pow>0.99){
						
							break;
					}
					fraction_high+=lambda2.get(i)/sum_lambda;
					lambda_high.add(lambda2.get(i));
				}
			}
		
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
			
			double[] cum_all=new double[lambda2.size()];
			pp=0;
			for (int i=0;i<cum_all.length;i++){
				cum_all[i]=pp;
				pp+=lambda2.get(i);
				
			}
			for (int i=0;i<cum_all.length;i++){
				cum_all[i]/=pp;
			}
		
			ArrayList<Double> s1=new ArrayList<Double>();
			ArrayList<Double> s3=new ArrayList<Double>();
			
			
			BinomialDistribution binom=new BinomialDistribution(NN,fraction_high);
			//System.out.println("SART");
			for (int k1=0;k1<10000;k1++){//10000
				int N_high=binom.sample();
				int[] meta_random=random_collision(N_high,cum_high);
				meta_random[0]+=NN-N_high;
				double sum1=0;
				for (int i=1;i<meta_random.length;i++){
					if(meta_random[i]>0){
						sum1+=Gamma.logGamma(i+1+1)*meta_random[i];
					}
				}
				s1.add(-sum1);
				//for (int i=0;i<meta_random.length;i++){
				//	System.out.print("	"+meta_random[i]);
				//}
				//System.out.println();
			}
			//System.out.println("END");
			
		
			
			
			for (int k2=0;k2<100000;k2++){//100000
				
				s3.add(offset+random_draw(NN,lambda_dist,err,lambda2));
			//	System.out.println(s3.get(s3.size()-1));
			}
			//System.out.println("END2");
			Collections.sort(s1);
			Collections.sort(s3);
			//System.out.println("END3");
			//System.out.println("Compare");
			/*double sum1=0;
			double sum2=0;
			int nn=0;
			for (int i=0;i<s1.size();i++){
				for (int j=0;j<s3.size();j++){
					if(s1.get(i)+s3.get(j)<=T){
						sum1+=Math.exp(s1.get(i)+s3.get(j));
						nn++;
					}
					sum2+=Math.exp(s1.get(i)+s3.get(j));
				///	System.out.println(i+"	"+j+"	"+Math.exp(s1.get(i)+s3.get(j)));
				}
			}*/
			
			
				int nn=0;
				for (int i=0;i<s3.size();i++){
					if(s3.get(i)+s1.get(0)>T){
						break;
					}
					for (int j=0;j<s1.size();j++){
						if(s3.get(i)+s1.get(j)>T){
							break;
						}
						nn++;
					}
				}
				return (double)(nn)/(double)(s1.size()*s3.size());
			
			
			
			//System.out.println("frac: "+(double)(nn)/(double)(s1.size()*s3.size()));
			
			
			//System.out.println("sum1: "+sum1);
			//System.out.println("sum2: "+sum2);
			//System.out.println("frac_sum: "+(sum1/sum2));
		/*System.out.println("SUM1");
			for (int i=0;i<s1.size();i++){
				System.out.println(Math.exp(s1.get(i)));
			}
			System.out.println("SUM2");
			for (int i=0;i<s3.size();i++){
				System.out.println(Math.exp(s3.get(i)));
			}
			*/
			/*
			System.out.println("Experimental");
			double[] cum=new double[lambda2.size()];
			double p=0;
			for (int i=0;i<cum.length;i++){
				cum[i]=p;
				p+=lambda2.get(i);
			}
			for (int i=0;i<cum.length;i++){
				cum[i]/=p;
			}
			for (int k=0;k<100;k++){
				int[] count=new int[cum.length];
				for (int i=0;i<NN;i++){
					double r=Math.random();
					for (int j=0;j<cum.length;j++){
						if(j<cum.length-1){
							if(cum[j]<=r&&r<cum[j+1]){
								count[j]++;
								break;
							}
						}
						else{
							if(cum[cum.length-1]<=r){
								count[cum.length-1]++;
								break;
							}
						}

					}
				}
				double ss=0;
				for (int i=0;i<count.length;i++){
					if(count[i]>0){
						ss+=Math.log(lambda2.get(i))*count[i];
					}
					if(count[i]>1){
						ss-=Gamma.logGamma(1+count[i]);
					}
					
				}
				System.out.println(Math.exp(ss));
			}
			*/
			
			
			
			
			
			
			
			
			
			
			/*
			double sum1=0;
			double sum2=0;
			double sum3=0;
			for (int i=0;i<s1.size();i++){
				sum1+=(s1.get(i));
			}
			for (int i=0;i<s3.size();i++){
				sum2+=(s3.get(i));
			}
			for (int i=0;i<s3.size();i++){
				if(s3.get(i)+s1.get(0)>T){
					break;
				}
				for (int j=0;j<s1.size();j++){
					if(s3.get(i)+s1.get(j)>T){
						break;
					}
				//	System.out.println(s3.get(i)+s1.get(j));
					sum3+=(s3.get(i)+s1.get(j));
				}
			}*/
			
			/*
			System.out.println("Sum1: "+sum1);
			System.out.println("Sum2: "+sum2);
			System.out.println("Sum3: "+sum3);
			
			return sum3-(s3.size()*sum1+s1.size()*sum2);
			*/
		}
		
		public double sign_quick(ArrayList<Double> lambda2, ArrayList<Integer> count2, int err){
			
			
			int NN=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
			}
			if(lambda2.size()==0||NN==0){
				return 1;
			}
			
			double T=0;
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
					//System.out.println(Math.log(lambda2.get(i)));
				}
				if(count2.get(i)>0){
					T-=Gamma.logGamma(count2.get(i)+1);
				}
			}
			double offset=NN*Math.log((double)(1)/(double)(lambda2.size()));
			T+=offset;
			
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			//int NN=500;
			
			Collections.sort(lambda2,comp);
			double sum_lambda=0;
			for (int i=0;i<lambda2.size();i++){
				sum_lambda+=lambda2.get(i);
			}
			
			double sum_pow=0;
			{
				int k=2;
				for (int i=0;i<lambda2.size();i++){
					sum_pow+=Math.pow(lambda2.get(i),k);
				}
			}
		
			

			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			{
				int k=2;
				double sum_local=0;
				for (int i=0;i<lambda2.size();i++){
					sum_local+=Math.pow(lambda2.get(i),k);
					if(sum_local/sum_pow>0.99){
							break;
					}
					fraction_high+=lambda2.get(i)/sum_lambda;
					lambda_high.add(lambda2.get(i));
				}
			}
		
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
			
			double[] cum_all=new double[lambda2.size()];
			pp=0;
			for (int i=0;i<cum_all.length;i++){
				cum_all[i]=pp;
				pp+=lambda2.get(i);
				
			}
			for (int i=0;i<cum_all.length;i++){
				cum_all[i]/=pp;
			}
		
			ArrayList<Double> s1=new ArrayList<Double>();
			ArrayList<Double> s3=new ArrayList<Double>();
			
			
			BinomialDistribution binom=new BinomialDistribution(NN,fraction_high);
			//System.out.println("SART");
			for (int k1=0;k1<100;k1++){//
				int N_high=binom.sample();
				int[] meta_random=random_collision(N_high,cum_high);
				meta_random[0]+=NN-N_high;
				double sum1=0;
				for (int i=1;i<meta_random.length;i++){
					if(meta_random[i]>0){
						sum1+=Gamma.logGamma(i+1+1)*meta_random[i];
					}
				}
				s1.add(-sum1);
				//for (int i=0;i<meta_random.length;i++){
				//	System.out.print("	"+meta_random[i]);
				//}
				//System.out.println();
			}
			//System.out.println("END");
			for (int k2=0;k2<1000;k2++){//
				
				s3.add(offset+random_draw(NN,lambda_dist,err,lambda2));
			}
			
			Collections.sort(s1);
			Collections.sort(s3);
			
			int n1=0;
			int n2=0;
			for (int k=0;k<1000;k++){
				if(s3.get((int)(s3.size()*random_conc.nextDouble()))+s1.get((int)(s1.size()*random_conc.nextDouble()))<=T){
					n1++;
				}
				n2++;
			}
			
			return (double)(n1)/(double)(n2);
				
			
		}
		
		
		
		public double sign_quick(ArrayList<Double> lambda2, ArrayList<Integer> count2, ArrayList<Double> probability){//int err, 
			
			if(probability.size()==0){
				return 1;
			}
			int NN=0;
			double sum_lambda=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
				sum_lambda+=lambda2.get(i);
			}
			if(lambda2.size()==0||NN==0){
				return 1;
			}
			if(probability.size()<=NN){
				return 0;
			}
			if(probability.get(NN)==0){
				return 0;
			}
			
			double T=Gamma.logGamma(NN+1);
			//double T1=0;
			//double T2=0;
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
					//T1+=count2.get(i)*Math.log(lambda2.get(i));
				}
				if(count2.get(i)>1){
					T-=Gamma.logGamma(count2.get(i)+1);
					//T2+=Gamma.logGamma(count2.get(i)+1);
				}
			}
			T-=Math.log(sum_lambda)*NN;
		//	T+=Math.log(probability.get(NN));
			
			
			
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			
			Collections.sort(lambda2,comp);
			
		
			double sum_sq=0;
			for (int i=0;i<lambda2.size();i++){
				sum_sq+=Math.pow(lambda2.get(i),2);
			}
		
			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),2);
				if(sum_local/sum_sq>0.99){
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
			}
			
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
			/*
			double[] cum_all=new double[lambda2.size()];
			pp=0;
			for (int i=0;i<cum_all.length;i++){
				cum_all[i]=pp;
				pp+=lambda2.get(i);
				
			}
			for (int i=0;i<cum_all.length;i++){
				cum_all[i]/=pp;
			}
			*/
			
			double[][] s1=new double[probability.size()-1][100]; 
			double[][] s3=new double[probability.size()-1][1000]; 
			
			int[][] count=new int[s1[0].length][lambda_high.size()];
			for (int i=0;i<s1[0].length;i++){
				for (int j=0;j<s1.length;j++){
					if(random_conc.nextDouble()>fraction_high){
						if(j==0){
							s1[j][i]=0;
						}
						else{
							s1[j][i]=s1[j-1][i];
						}
					}
					else{
						int r=random_draw(cum_high);
						if(count[i][r]==0){
							if(j==0){
								s1[j][i]=0;
							}
							else{
								s1[j][i]=s1[j-1][i];
							}
						}
						else if(count[i][r]==1){
							s1[j][i]=s1[j-1][i]+Gamma.logGamma(2+1);
						}
						else{
							s1[j][i]=s1[j-1][i]+Gamma.logGamma(count[i][r]+1+1)-Gamma.logGamma(count[i][r]+1);
						}
						count[i][r]++;
					}
				}
			}
			
			for (int i=0;i<s3[0].length;i++){
				for (int j=0;j<s3.length;j++){
					if(j==0){
						s3[j][i]=random_draw_log(lambda_dist);
					}
					else{
						s3[j][i]=s3[j-1][i]+random_draw_log(lambda_dist);
					}
				}
			}
			/*
			for (int i=0;i<s1.length;i++){
				System.out.print(""+i);
				for (int j=0;j<s1[i].length;j++){
					System.out.print("	"+s1[i][j]);
				}
				System.out.println();
			}
			System.exit(0);*/
			for (int i=0;i<s1.length;i++){
				Arrays.sort(s1[i]);
			}
			for (int i=0;i<s3.length;i++){
				Arrays.sort(s3[i]);
			}
			
			
			
			double ppp=0;
			double [] cum_probability=new double[probability.size()];
			double sum_prob=0;
			for (int i=0;i<probability.size();i++){
				sum_prob+=probability.get(i);
			}
			for (int i=0;i<probability.size();i++){
				cum_probability[i]=ppp/sum_prob;
				ppp+=probability.get(i);
			}
			
			
//			for (int i=0;i<probability.size();i++){
//				System.out.println(i+"	"+probability.get(i));
//			}
			
			int n1=0;
			int n2=0;
//			System.out.println("T	"+T);
//			System.out.println("T	"+NN);
			for (int k=0;k<1000;k++){
				
				int M=random_draw(cum_probability);
				//int M=NN;
		//		System.out.println("exp	"+M);
				double exp=Gamma.logGamma(M+1)-M*Math.log(sum_lambda);//Math.log(probability.get(M))+
				if(M>0){
					exp+=-s1[M-1][(int)(random_conc.nextDouble()*s1[0].length)]+s3[M-1][(int)(random_conc.nextDouble()*s3[0].length)];
				}
				if(exp<=T){
					n1++;
				}
				n2++;
	//			System.out.println("exp	"+exp);
		//		System.out.println(T+"	"+exp);
				//System.out.println(s1[M-1][(int)(random_conc.nextDouble()*s1[0].length)]+"	"+T2);
			}
	//		System.exit(0);
			return (double)(n1)/(double)(n2);
			
			
			/*
			int n1=0;
			int n2=0;
			for (int k=0;k<1000;k++){
				if(s3.get((int)(s3.size()*random_conc.nextDouble()))+s1.get((int)(s1.size()*random_conc.nextDouble()))<=T){
					n1++;
				}
				n2++;
			}
			
			return (double)(n1)/(double)(n2);
			*/	
			
		}
		public double sign_quickREVISED(ArrayList<Double> lambda2, ArrayList<Integer> count2, ArrayList<Double> probability){//int err, 
			//System.out.println("STARTING Sign Revised");
			if(probability.size()==0){
				return 1;
			}
			int NN=0;
			double sum_lambda=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
				sum_lambda+=lambda2.get(i);
			}
			//System.out.println("NN: "+NN);
			if(lambda2.size()==0||NN==0){
				return 1;
			}
			if(probability.size()<=NN){
				return 0;
			}
			if(probability.get(NN)==0){
				return 0;
			}
			
			double T=Gamma.logGamma(NN+count2.size()+1);
		
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
					//T1+=count2.get(i)*Math.log(lambda2.get(i));
				}
				if(count2.get(i)>1){
					T-=Gamma.logGamma(count2.get(i)+1);
					//T2+=Gamma.logGamma(count2.get(i)+1);
				}
			}
			T-=Math.log(sum_lambda)*NN;
			
			
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			
			Collections.sort(lambda2,comp);
			
		
			double sum_sq=0;
			for (int i=0;i<lambda2.size();i++){
				sum_sq+=Math.pow(lambda2.get(i),2);
			}
		
			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),2);
				if(sum_local/sum_sq>0.99){
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
			}
			
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
		
			
			double[][] s1=new double[probability.size()-1][100]; 
			double[][] s3=new double[probability.size()-1][1000]; 
			
			int[][] count=new int[s1[0].length][lambda_high.size()];
			for (int i=0;i<s1[0].length;i++){
				for (int j=0;j<s1.length;j++){
					if(random_conc.nextDouble()>fraction_high){
						if(j==0){
							s1[j][i]=0;
						}
						else{
							s1[j][i]=s1[j-1][i];
						}
					}
					else{
						int r=random_draw(cum_high);
						if(count[i][r]==0){
							if(j==0){
								s1[j][i]=0;
							}
							else{
								s1[j][i]=s1[j-1][i];
							}
						}
						else if(count[i][r]==1){
							s1[j][i]=s1[j-1][i]+Gamma.logGamma(2+1);
						}
						else{
							s1[j][i]=s1[j-1][i]+Gamma.logGamma(count[i][r]+1+1)-Gamma.logGamma(count[i][r]+1);
						}
						count[i][r]++;
					}
				}
			}
			
			for (int i=0;i<s3[0].length;i++){
				for (int j=0;j<s3.length;j++){
					if(j==0){
						s3[j][i]=random_draw_log(lambda_dist);
					}
					else{
						s3[j][i]=s3[j-1][i]+random_draw_log(lambda_dist);
					}
				}
			}
			
			for (int i=0;i<s1.length;i++){
				Arrays.sort(s1[i]);
			}
			for (int i=0;i<s3.length;i++){
				Arrays.sort(s3[i]);
			}
			
			
			
			double ppp=0;
			double [] cum_probability=new double[probability.size()];
			double sum_prob=0;
			for (int i=0;i<probability.size();i++){
				sum_prob+=probability.get(i);
			}
			for (int i=0;i<probability.size();i++){
				cum_probability[i]=ppp/sum_prob;
				ppp+=probability.get(i);
			}
			
			

			int n1=0;
			int n2=0;
			
			double frac_part=0;
			for (int i=0;i<Math.min(NN,probability.size());i++){
				frac_part+=probability.get(i);
			}

			for (int k=0;k<1000;k++){
				
				int M=random_draw(cum_probability,NN);
				
				double exp=Gamma.logGamma(M+count2.size()+1)-M*Math.log(sum_lambda);//Math.log(probability.get(M))+
				if(M>0){
					exp+=-s1[M-1][(int)(random_conc.nextDouble()*s1[0].length)]+s3[M-1][(int)(random_conc.nextDouble()*s3[0].length)];
				}
				if(exp<=T){
					n1++;
				}
				n2++;
			}
			return (double)(n1)/(double)(n2)*(1-frac_part);
			
	
		}
		
		
		
		public double sign_quickREVISED(int count_ext, ArrayList<Double> lambda2, ArrayList<Integer> count2, ArrayList<Double> probability, int err){//int err, 
			//TODO:err
			if(lambda2.size()==0||count_ext==0||probability.size()==0){
				return 1;
			}
			if(count_ext>=probability.size()||probability.get(count_ext)==0){
				return 0;
			}
			//System.out.println("STARTING Sign Revised");
			if(probability.size()==0){
				return 1;
			}
			int NN=0;
			double sum_lambda=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
				sum_lambda+=lambda2.get(i);
			}
			//System.out.println("NN: "+NN);
			
			if(probability.size()<=NN){
				return 0;
			}
			if(probability.get(NN)==0){
				return 0;
			}
			
			double T=Gamma.logGamma(NN+count2.size()+1);
		
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
					//T1+=count2.get(i)*Math.log(lambda2.get(i));
				}
				if(count2.get(i)>1){
					T-=Gamma.logGamma(count2.get(i)+1);
					//T2+=Gamma.logGamma(count2.get(i)+1);
				}
			}
			T-=Math.log(sum_lambda)*NN;
			T+=Math.log(probability.get(NN));
			
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			
			Collections.sort(lambda2,comp);
			
		
			double sum_sq=0;
			for (int i=0;i<lambda2.size();i++){
				sum_sq+=Math.pow(lambda2.get(i),2);
			}
		
			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),2);
				if(sum_local/sum_sq>0.99){
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
				if(lambda_high.size()>50000){
					return 1;
				}
			}
			
			
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
		
			
			double[][] s1=new double[probability.size()-1][100]; 
			double[][] s3=new double[probability.size()-1][1000]; 
			
			int[][] count=new int[s1[0].length][lambda_high.size()];
			for (int i=0;i<s1[0].length;i++){
				for (int j=0;j<s1.length;j++){
					if(random_conc.nextDouble()>fraction_high){
						if(j==0){
							s1[j][i]=0;
						}
						else{
							s1[j][i]=s1[j-1][i];
						}
					}
					else{
						int r=random_draw(cum_high);
						if(count[i][r]==0){
							if(j==0){
								s1[j][i]=0;
							}
							else{
								s1[j][i]=s1[j-1][i];
							}
						}
						else if(count[i][r]==1){
							s1[j][i]=s1[j-1][i]+Gamma.logGamma(2+1);
						}
						else{
							s1[j][i]=s1[j-1][i]+Gamma.logGamma(count[i][r]+1+1)-Gamma.logGamma(count[i][r]+1);
						}
						count[i][r]++;
					}
				}
			}
			
			for (int i=0;i<s3[0].length;i++){
				for (int j=0;j<s3.length;j++){
					if(j==0){
						s3[j][i]=random_draw_log(lambda_dist);
					}
					else{
						s3[j][i]=s3[j-1][i]+random_draw_log(lambda_dist);
					}
				}
			}
			
			for (int i=0;i<s3.length;i++){
				for (int j=0;j<s3[i].length;j++){
					for (int k=0;k<err;k++){
						s3[i][j]+=Math.log(lambda2.get((int)(random_conc.nextDouble()*lambda2.size())));
					}
				}
			}
			
			for (int i=0;i<s1.length;i++){
				Arrays.sort(s1[i]);
			}
			for (int i=0;i<s3.length;i++){
				Arrays.sort(s3[i]);
			}
			
			
			
			double ppp=0;
			double [] cum_probability=new double[probability.size()];
			double sum_prob=0;
			for (int i=0;i<probability.size();i++){
				sum_prob+=probability.get(i);
			}
			for (int i=0;i<probability.size();i++){
				cum_probability[i]=ppp/sum_prob;
				ppp+=probability.get(i);
			}
			
			

			int n1=0;
			int n2=0;
			
			double frac_part=0;
			for (int i=0;i<Math.min(count_ext,probability.size());i++){
				frac_part+=probability.get(i);
			}

			for (int k=0;k<1000;k++){
				
				int M=random_draw(cum_probability,count_ext);
				
				double exp=Gamma.logGamma(M+count2.size()+1)-M*Math.log(sum_lambda)+Math.log(probability.get(M));
				if(M>0){
					exp+=-s1[M-1][(int)(random_conc.nextDouble()*s1[0].length)]+s3[M-1][(int)(random_conc.nextDouble()*s3[0].length)];
				}
				if(exp<=T){
					n1++;
				}
				n2++;
			}
			return (double)(n1)/(double)(n2)*(1-frac_part);
			
	
		}
		
		
		public double sign_long(ArrayList<Double> lambda2, ArrayList<Integer> count2, ArrayList<Double> probability){//int err, 
//			System.out.println("Entering into long");
			if(probability.size()==0){
				return 1;
			}
			int NN=0;
			double sum_lambda=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
				sum_lambda+=lambda2.get(i);
			}
			if(lambda2.size()==0||NN==0){
				return 1;
			}
			if(probability.size()<=NN){
				System.out.println("Probability too short");
				System.exit(0);
				return 0;
			}
			if(probability.get(NN)==0){
				return 0;
			}
			
			double T=Gamma.logGamma(NN+1);
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
				}
				if(count2.get(i)>1){
					T-=Gamma.logGamma(count2.get(i)+1);
				}
			}
			T-=Math.log(sum_lambda)*NN;
			//T+=Math.log(probability.get(NN));
			
			
			
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			
			Collections.sort(lambda2,comp);
			
		
			double sum_sq=0;
			for (int i=0;i<lambda2.size();i++){
				sum_sq+=Math.pow(lambda2.get(i),2);
			}
		
			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),2);
				if(sum_local/sum_sq>0.99){
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
			}
			
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
//			System.out.println("generating arrays");
			
			double[][] s1=new double[probability.size()-1][10000]; 
			double[][] s3=new double[probability.size()-1][100000]; 
			
//			System.out.println("generating arrays done");
			
			
			int[][] count=new int[s1[0].length][lambda_high.size()];
			for (int i=0;i<s1[0].length;i++){
				for (int j=0;j<s1.length;j++){
					if(random_conc.nextDouble()>fraction_high){
						if(j==0){
							s1[j][i]=0;
						}
						else{
							s1[j][i]=s1[j-1][i];
						}
					}
					else{
						int r=random_draw(cum_high);
						if(count[i][r]==0){
							if(j==0){
								s1[j][i]=0;
							}
							else{
								s1[j][i]=s1[j-1][i];
							}
						}
						else if(count[i][r]==1){
							s1[j][i]=s1[j-1][i]-Gamma.logGamma(2+1);
						}
						else{
							s1[j][i]=s1[j-1][i]-Gamma.logGamma(count[i][r]+1+1)+Gamma.logGamma(count[i][r]+1);
						}
						count[i][r]++;
					}
				}
			}
//			System.out.println("s1 done");
			
			for (int i=0;i<s3[0].length;i++){
				for (int j=0;j<s3.length;j++){
					if(j==0){
						s3[j][i]=random_draw_log(lambda_dist);
					}
					else{
						s3[j][i]=s3[j-1][i]+random_draw_log(lambda_dist);
					}
				}
			}
//			System.out.println("s3 done");

			for (int i=0;i<s1.length;i++){
				Arrays.sort(s1[i]);
			}
			for (int i=0;i<s3.length;i++){
				Arrays.sort(s3[i]);
			}
			
			double sign=0;
			for (int M=0;M<probability.size();M++){
				if(probability.get(M)==0){
					continue;
				}
				double start=+Gamma.logGamma(M+1)-M*Math.log(sum_lambda);//Math.log(probability.get(M))
				double frac=1;
				if(M>0){
					int nn=0;
					for (int i=0;i<1000;i++){
						if(s3[M-1][(int)(random_conc.nextDouble()*s3[M-1].length)]+s1[M-1][(int)(s1[M-1].length*random_conc.nextDouble())]+start<=T){
							nn++;
						}
					}
					double f=(double)(nn)/1000.0;
					if(f<=0.05){
						//System.out.println("Intense "+f);
						nn=0;
						for (int i=0;i<s3[M-1].length;i++){
							if(s3[M-1][i]+s1[M-1][0]+start>T){
								break;
							}
							for (int j=0;j<s1[M-1].length;j++){
								if(s3[M-1][i]+s1[M-1][j]+start>T){
									break;
								}
								nn++;
							}
						}
						frac=(double)(nn)/(double)(s1[M-1].length*s3[M-1].length);
					}
					else{
						frac=f;
					}
					
				}
				else{
					if(Math.log(probability.get(M))<=T){
						frac=1;
					}
					else{
						frac=0;
					}
				}
//				System.out.println(M+"	"+NN+"	"+frac+"	"+probability.get(M));
				
				//System.out.println(frac+"");
				sign+=probability.get(M)*frac;
			}
//			System.out.println(sign);
			
			return sign;
/*
			
			int nn=0;
			
			
			for (int i=0;i<s3.size();i++){
				if(s3.get(i)+s1.get(0)>T){
					break;
				}
				for (int j=0;j<s1.size();j++){
					if(s3.get(i)+s1.get(j)>T){
						break;
					}
					nn++;
				}
			}
			return (double)(nn)/(double)(s1.size()*s3.size());
		
		
			*/
			/*
			double ppp=0;
			double [] cum_probability=new double[probability.size()];
			double sum_prob=0;
			for (int i=0;i<probability.size();i++){
				sum_prob+=probability.get(i);
			}
			for (int i=0;i<probability.size();i++){
				cum_probability[i]=ppp/sum_prob;
				ppp+=probability.get(i);
			}

			
			int n1=0;
			int n2=0;

			for (int k=0;k<1000;k++){
				
				int M=random_draw(cum_probability);
				double exp=Math.log(probability.get(M))+Gamma.logGamma(M+1)-M*Math.log(sum_lambda);
				if(M>0){
					exp+=-s1[M-1][(int)(random_conc.nextDouble()*s1[0].length)]+s3[M-1][(int)(random_conc.nextDouble()*s3[0].length)];
				}
				if(exp<=T){
					n1++;
				}
				n2++;
			}
			return (double)(n1)/(double)(n2);*/
			
		}
		
		
		
		public double sign_longREVISED(ArrayList<Double> lambda2, ArrayList<Integer> count2, ArrayList<Double> probability){//int err, 
			if(probability.size()==0){
				return 1;
			}
			int NN=0;
			double sum_lambda=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
				sum_lambda+=lambda2.get(i);
			}
			if(lambda2.size()==0||NN==0){
				return 1;
			}
			if(probability.size()<=NN){
				System.out.println("Probability too short");
				System.exit(0);
				return 0;
			}
			if(probability.get(NN)==0){
				return 0;
			}
			
			double T=Gamma.logGamma(NN+count2.size()+1);
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
				}
				if(count2.get(i)>1){
					T-=Gamma.logGamma(count2.get(i)+1);
				}
			}
			T-=Math.log(sum_lambda)*NN;
		
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			
			Collections.sort(lambda2,comp);
			
		
			double sum_sq=0;
			for (int i=0;i<lambda2.size();i++){
				sum_sq+=Math.pow(lambda2.get(i),2);
			}
		
			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),2);
				if(sum_local/sum_sq>0.99){
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
			}
			
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
//			System.out.println("generating arrays");
			
			double[][] s1=new double[probability.size()-1][10000]; 
			double[][] s3=new double[probability.size()-1][100000]; 
			
//			System.out.println("generating arrays done");
			
			
			int[][] count=new int[s1[0].length][lambda_high.size()];
			for (int i=0;i<s1[0].length;i++){
				for (int j=0;j<s1.length;j++){
					if(random_conc.nextDouble()>fraction_high){
						if(j==0){
							s1[j][i]=0;
						}
						else{
							s1[j][i]=s1[j-1][i];
						}
					}
					else{
						int r=random_draw(cum_high);
						if(count[i][r]==0){
							if(j==0){
								s1[j][i]=0;
							}
							else{
								s1[j][i]=s1[j-1][i];
							}
						}
						else if(count[i][r]==1){
							s1[j][i]=s1[j-1][i]-Gamma.logGamma(2+1);
						}
						else{
							s1[j][i]=s1[j-1][i]-Gamma.logGamma(count[i][r]+1+1)+Gamma.logGamma(count[i][r]+1);
						}
						count[i][r]++;
					}
				}
			}
//			System.out.println("s1 done");
			
			for (int i=0;i<s3[0].length;i++){
				for (int j=0;j<s3.length;j++){
					if(j==0){
						s3[j][i]=random_draw_log(lambda_dist);
					}
					else{
						s3[j][i]=s3[j-1][i]+random_draw_log(lambda_dist);
					}
				}
			}
//			System.out.println("s3 done");

			for (int i=0;i<s1.length;i++){
				Arrays.sort(s1[i]);
			}
			for (int i=0;i<s3.length;i++){
				Arrays.sort(s3[i]);
			}
			
			double sign=0;
			for (int M=NN;M<probability.size();M++){
				if(probability.get(M)==0){
					continue;
				}
				double start=+Gamma.logGamma(M+count2.size()+1)-M*Math.log(sum_lambda);//Math.log(probability.get(M))
				double frac=1;
				if(M>0){
					int nn=0;
					for (int i=0;i<1000;i++){
						if(s3[M-1][(int)(random_conc.nextDouble()*s3[M-1].length)]+s1[M-1][(int)(s1[M-1].length*random_conc.nextDouble())]+start<=T){
							nn++;
						}
					}
					double f=(double)(nn)/1000.0;
					if(f<=0.05){
						//System.out.println("Intense "+f);
						nn=0;
						for (int i=0;i<s3[M-1].length;i++){
							if(s3[M-1][i]+s1[M-1][0]+start>T){
								break;
							}
							for (int j=0;j<s1[M-1].length;j++){
								if(s3[M-1][i]+s1[M-1][j]+start>T){
									break;
								}
								nn++;
							}
						}
						frac=(double)(nn)/(double)(s1[M-1].length*s3[M-1].length);
					}
					else{
						frac=f;
					}
					
				}
				else{
					if(Math.log(probability.get(M))<=T){
						frac=1;
					}
					else{
						frac=0;
					}
				}
				sign+=probability.get(M)*frac;
				
			}
			return sign;

			
		}
		
		
		
		public double sign_longREVISED(int count_ext, ArrayList<Double> lambda2, ArrayList<Integer> count2, ArrayList<Double> probability, int err){//int err, 
			//TODO err
			if(lambda2.size()==0||count_ext==0||probability.size()==0){
				return 1;
			}
			if(count_ext>=probability.size()||probability.get(count_ext)==0){
				return 0;
			}
			
			if(probability.size()==0){
				return 1;
			}
			int NN=0;
			double sum_lambda=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
				sum_lambda+=lambda2.get(i);
			}
			if(lambda2.size()==0||NN==0){
				return 1;
			}
			if(probability.size()<=NN){
				System.out.println("Probability too short");
				System.exit(0);
				return 0;
			}
			if(probability.get(NN)==0){
				return 0;
			}
			
			double T=Gamma.logGamma(NN+count2.size()+1);
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
				}
				if(count2.get(i)>1){
					T-=Gamma.logGamma(count2.get(i)+1);
				}
			}
			T-=Math.log(sum_lambda)*NN;
			T+=Math.log(probability.get(NN));
			
			
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			
			Collections.sort(lambda2,comp);
			
		
			double sum_sq=0;
			for (int i=0;i<lambda2.size();i++){
				sum_sq+=Math.pow(lambda2.get(i),2);
			}
		
			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),2);
				if(sum_local/sum_sq>0.99){
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
				if(lambda_high.size()>50000){
					return 1;
				}
			}
			
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
//			System.out.println("generating arrays");
			
			double[][] s1=new double[probability.size()-1][10000]; 
			double[][] s3=new double[probability.size()-1][100000]; 
			
//			System.out.println("generating arrays done");
			
			//System.out.println(lambda_high.size());//TODO: 
			int[][] count=new int[s1[0].length][lambda_high.size()];
			for (int i=0;i<s1[0].length;i++){
				for (int j=0;j<s1.length;j++){
					if(random_conc.nextDouble()>fraction_high){
						if(j==0){
							s1[j][i]=0;
						}
						else{
							s1[j][i]=s1[j-1][i];
						}
					}
					else{
						int r=random_draw(cum_high);
						if(count[i][r]==0){
							if(j==0){
								s1[j][i]=0;
							}
							else{
								s1[j][i]=s1[j-1][i];
							}
						}
						else if(count[i][r]==1){
							s1[j][i]=s1[j-1][i]-Gamma.logGamma(2+1);
						}
						else{
							s1[j][i]=s1[j-1][i]-Gamma.logGamma(count[i][r]+1+1)+Gamma.logGamma(count[i][r]+1);
						}
						count[i][r]++;
					}
				}
			}
//			System.out.println("s1 done");
			
			for (int i=0;i<s3[0].length;i++){
				for (int j=0;j<s3.length;j++){
					if(j==0){
						s3[j][i]=random_draw_log(lambda_dist);
					}
					else{
						s3[j][i]=s3[j-1][i]+random_draw_log(lambda_dist);
					}
				}
			}
//			System.out.println("s3 done");
			for (int i=0;i<s3.length;i++){
				for (int j=0;j<s3[i].length;j++){
					for (int k=0;k<err;k++){
						s3[i][j]+=Math.log(lambda2.get((int)(random_conc.nextDouble()*lambda2.size())));
					}
				}
			}
			
			
			for (int i=0;i<s1.length;i++){
				Arrays.sort(s1[i]);
			}
			for (int i=0;i<s3.length;i++){
				Arrays.sort(s3[i]);
			}
			
			double sign=0;
			for (int M=count_ext;M<probability.size();M++){
				if(probability.get(M)==0){
					continue;
				}
				double start=+Gamma.logGamma(M+count2.size()+1)-M*Math.log(sum_lambda)+Math.log(probability.get(M));
				double frac=1;
				if(M>0){
					int nn=0;
					for (int i=0;i<1000;i++){
						if(s3[M-1][(int)(random_conc.nextDouble()*s3[M-1].length)]+s1[M-1][(int)(s1[M-1].length*random_conc.nextDouble())]+start<=T){
							nn++;
						}
					}
					double f=(double)(nn)/1000.0;
					if(f<=0.05){
						//System.out.println("Intense "+f);
						nn=0;
						for (int i=0;i<s3[M-1].length;i++){
							if(s3[M-1][i]+s1[M-1][0]+start>T){
								break;
							}
							for (int j=0;j<s1[M-1].length;j++){
								if(s3[M-1][i]+s1[M-1][j]+start>T){
									break;
								}
								nn++;
							}
						}
						frac=(double)(nn)/(double)(s1[M-1].length*s3[M-1].length);
					}
					else{
						frac=f;
					}
					
				}
				else{
					if(Math.log(probability.get(M))<=T){
						frac=1;
					}
					else{
						frac=0;
					}
				}
				sign+=probability.get(M)*frac;
				
			}
			return sign;

			
		}
		
		
		/*
		public double sign(ArrayList<Double> lambda2, ArrayList<Integer> count2 ,double rate, double cov){
			double sign_quick=sign_quick(lambda2,count2, rate,cov);
			if(sign_quick<0.05){
				return sign_long(lambda2,count2, rate,cov);
			}
			else{
				return sign_quick;
			}
		}*/
		
		public int random_draw(double[] cum){
			double r=random_conc.nextDouble();
			for (int i=0;i<cum.length;i++){
				if(i<cum.length-1){
					if(cum[i]<=r&&r<cum[i+1]){
						return i;
					}
				}
				else{
					if(cum[i]<=r){
						return i;
					}
				}
			}
			return -1;
		}
		
		
		public int random_draw(double[] cum, int min){
			
			double r=(1-cum[min])*random_conc.nextDouble()+cum[min];
			for (int i=0;i<cum.length;i++){
				if(i<cum.length-1){
					if(cum[i]<=r&&r<cum[i+1]){
						return i;
					}
				}
				else{
					if(cum[i]<=r){
						return i;
					}
				}
			}
			return -1;
		}
				
		
		public double random_draw(int N, double[] lambda, int N_err, ArrayList<Double> lambda2){
			if(N==0){
				return 0;
			}
			double sum=0;
			for (int i=0;i<N;i++){
				sum+=lambda[(int)(lambda.length*random_conc.nextDouble())];
			}
			for (int i=0;i<N_err;i++){
				sum+=Math.log(lambda2.get((int)(random_conc.nextDouble()*lambda2.size())));
			}
			return sum;
		}
		
		public double hotspot_sign(ArrayList<double[]> pos){
			int max_count=0;
			double lambda_sum=0;
			int count_sum=0;
			
			for (int i=0;i<pos.size();i++){
				if(pos.get(i)[0]>max_count){
					max_count=(int)(pos.get(i)[0]);
				}
				count_sum+=(int)(pos.get(i)[0]);
				lambda_sum+=pos.get(i)[1];
			}
			
			if(max_count>=4&&(double)(max_count)/(double)(count_sum)>=0.05){
				double p=1;
				for (int i=0;i<pos.size();i++){
					p*=new BinomialDistribution(count_sum,pos.get(i)[1]/lambda_sum).cumulativeProbability(max_count-1);
				}
				return 1-p;
			}
			else{
				return 1;
			}
		}
		
		
		public double random_draw_log(double[] lambda){
			
			return lambda[(int)(lambda.length*random_conc.nextDouble())];
			
			
		}
		
		/*
		public int[] random(int N, ArrayList<Double> lambda){
			double[] cum=new double[lambda.size()];
			double sum_lambda=0;
			for (int i=0;i<lambda.size();i++){
				cum[i]=sum_lambda;
				sum_lambda+=lambda.get(i);
			}
			int[] random=new int[lambda.size()];
			for (int k=0;k<N;k++){
				double r=random_conc.nextDouble();
				for (int j=0;j<cum.length;j++){
					if(j<cum.length-1){
						if(cum[j]/sum_lambda<=r&&r<cum[j+1]/sum_lambda){
							random[j]++;
							break;
						}
					}
					else{
						if(cum[cum.length-1]/sum_lambda<=r){
							random[j]++;
							break;
						}
					}
					
				}
			}
			return random;
		}*/
		
		
		public void run(){
			System.out.println("start "+chr[c]);
			random_conc= ThreadLocalRandom.current();
			try{
				FileInputStream in=new FileInputStream(file_count+entity+"_Chr"+chr[c]+".txt");
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				
				FileInputStream in2=new FileInputStream(file_amino+chr[c]+".txt");
				DataInputStream inn2=new DataInputStream(in2);
				BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
				
				
				ArrayList<Integer> position=new ArrayList<Integer>();
				ArrayList<double[]> lambda=new ArrayList<double[]>();
				ArrayList<String> nucl=new ArrayList<String>();
				ArrayList<Double>coverage=new ArrayList<Double>();
				ArrayList<int[]> label=new ArrayList<int[]>();
				ArrayList<int[]> count=new ArrayList<int[]>();
				ArrayList<String> amino_acid=new ArrayList<String>();
				
				for (int nn=0;nn<genes[c].size();nn++){
					int ii=0;
					if(position.size()>0){
						while(ii<position.size()&&position.get(ii)<genes[c].get(nn).start){
							ii++;
						}
					}
					if(ii>0){
						for (int i=ii-1;i>=0;i--){
							nucl.remove(i);
							lambda.remove(i);
							position.remove(i);
							coverage.remove(i);
							count.remove(i);
							label.remove(i);
							amino_acid.remove(i);
						}
					}
					
					String s="";
					while ((s=input.readLine())!=null){
						String[] t=s.split("	");
						String[] t2=input2.readLine().split("	");
						if(!t[0].equals(t2[0])){
							System.out.println("chr "+chr[c]+": count and amino acid files non-compatible");
							System.exit(0);
						}
						
						if(t2.length>2){
							amino_acid.add(t2[2]+t2[3]);
						}
						else{
							amino_acid.add("-");
						}
						
						position.add(Integer.parseInt(t[0]));
						nucl.add(t[1]);
						coverage.add(Double.parseDouble(t[2]));
						String[] tt=t[3].split(";");
						label.add(new int[]{Integer.parseInt(tt[0]),Integer.parseInt(tt[1]),Integer.parseInt(tt[2])});
						tt=t[4].split(";");
						lambda.add(new double[]{Double.parseDouble(tt[0]),Double.parseDouble(tt[1]),Double.parseDouble(tt[2])});
						tt=t[5].split(";");
						count.add(new int[]{Integer.parseInt(tt[0]),Integer.parseInt(tt[1]),Integer.parseInt(tt[2])});
					
						if(Integer.parseInt(t[0])>=genes[c].get(nn).end){
							break;
						}
					}
					
					
					ArrayList<Double> lambda2=new ArrayList<Double>();
					ArrayList<Integer> count2=new ArrayList<Integer>();
					ArrayList<Double> coverage2=new ArrayList<Double>();
					ArrayList<Double> lambda2_syn=new ArrayList<Double>();
					ArrayList<Integer> count2_syn=new ArrayList<Integer>();
					ArrayList<Double> coverage2_syn=new ArrayList<Double>();
					ArrayList<double[]> pos=new ArrayList<double[]>();
					ArrayList<double[]> pos_syn=new ArrayList<double[]>();
					ArrayList<String> aa=new ArrayList<String>();
					ArrayList<String> aa_syn=new ArrayList<String>();
					
					for (int i=0;i<position.size();i++){
						if(!genes[c].get(nn).contains(position.get(i))){
							continue;
						}
						for (int j=0;j<3;j++){
							if(label.get(i)[j]==1&&coverage.get(i)>0.5){
								pos.add(new double[]{count.get(i)[j],lambda.get(i)[j]});
								aa.add(amino_acid.get(i));
							}
							else if(label.get(i)[j]==1&&coverage.get(i)>0.5){
								pos_syn.add(new double[]{count.get(i)[j],lambda.get(i)[j]});
								aa_syn.add(amino_acid.get(i));
							}
						}
						
						
						if(coverage.get(i)>=0.05){//0.3
							
							
							if(label.get(i)[0]!=0){
								//lambda2.add(coverage.get(i)*lambda.get(i)[0]);
								lambda2.add(lambda.get(i)[0]);
								
								count2.add(count.get(i)[0]);
								coverage2.add(coverage.get(i));
							}
							else{
								//lambda2_syn.add(coverage.get(i)*lambda.get(i)[0]);
								lambda2_syn.add(lambda.get(i)[0]);
								
								count2_syn.add(count.get(i)[0]);
								coverage2_syn.add(coverage.get(i));
							}
							if(label.get(i)[1]!=0){
								//lambda2.add(coverage.get(i)*lambda.get(i)[1]);
								lambda2.add(lambda.get(i)[1]);
								
								count2.add(count.get(i)[1]);
								coverage2.add(coverage.get(i));
							}
							else{
								//lambda2_syn.add(coverage.get(i)*lambda.get(i)[1]);
								lambda2_syn.add(lambda.get(i)[1]);
								
								count2_syn.add(count.get(i)[1]);
								coverage2_syn.add(coverage.get(i));
							}
							if(label.get(i)[2]!=0){
								//lambda2.add(coverage.get(i)*lambda.get(i)[2]);
								lambda2.add(lambda.get(i)[2]);
								
								count2.add(count.get(i)[2]);
								coverage2.add(coverage.get(i));
							}
							else{
								//lambda2_syn.add(coverage.get(i)*lambda.get(i)[2]);
								lambda2_syn.add(lambda.get(i)[2]);
								
								count2_syn.add(count.get(i)[2]);
								coverage2_syn.add(coverage.get(i));
							}
						}
						
	
					}
					
					
					
					
					
					ArrayList<double[]> x=sort(lambda2,count2);
					lambda2=new ArrayList<Double>();
					count2=new ArrayList<Integer>();
					for (int i=0;i<x.size();i++){
						lambda2.add(x.get(i)[0]);
						count2.add((int)(x.get(i)[1]));
					}
					x=sort(lambda2_syn,count2_syn);
					lambda2_syn=new ArrayList<Double>();
					count2_syn=new ArrayList<Integer>();
					for (int i=0;i<x.size();i++){
						lambda2_syn.add(x.get(i)[0]);
						count2_syn.add((int)(x.get(i)[1]));
					}
					
					int err=3;//TODO: 3 ??
					
					
					//if(genes[c].get(nn).name.equals("LOXHD1")){
					
					if(genes[c].get(nn).cov*genes[c].get(nn).cov_syn==0){
						genes[c].get(nn).sign_vector=1;
						genes[c].get(nn).sign_hotspot=1;//signREVISED(lambda2,count2,genes[c].get(nn).prob_nonsyn,err);//err,
						genes[c].get(nn).sign_combined=1;//err,
						genes[c].get(nn).sign_cbase=1;
						genes[c].get(nn).sign_vector_syn=1;
						genes[c].get(nn).sign_hotspot_syn=1;
					
					}
					else{
						genes[c].get(nn).sign_vector=sign(lambda2,count2,err);
						genes[c].get(nn).sign_hotspot=hotspot_sign(summarize(pos,aa));//signREVISED(lambda2,count2,genes[c].get(nn).prob_nonsyn,err);//err,
						genes[c].get(nn).sign_combined=signREVISED(genes[c].get(nn).count,lambda2,count2,genes[c].get(nn).prob_nonsyn,err);//err,
						genes[c].get(nn).sign_cbase=cum_prob(genes[c].get(nn).count,genes[c].get(nn).prob_nonsyn);
						genes[c].get(nn).sign_vector_syn=sign(lambda2_syn,count2_syn,err);
						genes[c].get(nn).sign_hotspot_syn=hotspot_sign(summarize(pos_syn,aa_syn));
									
					}
					
			
				}
				input.close();
				input2.close();
				
			}
			catch(Exception e){
				StackTraceElement[] aa=e.getStackTrace();
				for (int i=0;i<aa.length;i++){
					System.out.println(i+"	"+aa[i].getLineNumber());
				}
				System.out.println(e);
			}
			System.out.println("done "+chr[c]);
			done=true;
		}
		
		public ArrayList<double[]> summarize(ArrayList<double[]> pos, ArrayList<String> aa){
			boolean[] taken=new boolean[pos.size()];
			ArrayList<double[]> pos2=new ArrayList<double[]>();
			for (int i=0;i<pos.size();i++){
				if(taken[i]){
					continue;
				}
				if(aa.get(i).equals("-")){
					pos2.add(pos.get(i));
					taken[i]=true;
					continue;
				}
				
				double[] a=new double[2];
				for (int j=0;j<9;j++){
					if(i+j>=pos.size()){
						break;
					}
					if(!aa.get(i+j).equals("-")&&aa.get(i+j).equals(aa.get(i))){
						a[0]+=pos.get(i+j)[0];
						a[1]+=pos.get(i+j)[1];
						taken[i+j]=true;
					}
				}
				pos2.add(a);
			}
			return pos2;
		}
		
		
		public static double cum_prob(int count, ArrayList<Double> prob){
			if(prob.size()==0){
				return 1;
			}
			if(count==0){
				return 1;
			}
			double sum=0;
			for (int i=0;i<Math.min(prob.size(),count);i++){
				sum+=prob.get(i);
			}
			return 1-sum;
		}
	}
	
	/*
	public static int error(ArrayList<Double> lambda2, ArrayList<Integer> count2, int err){
		
		int n=count2.size()-1;
		int m=0;
		do{
			boolean found=false;
			for (;n>=0;n--){
				if(count2.get(n)>0){
					count2.set(n, 0);
					found=true;
					break;
				}
			}
			if(!found){
				break;
			}
			m++;
		}while(sign(lambda2,count2,err)<0.05);
		
		return m;
	}*/
	
	/*
	public static double error_sign(ArrayList<Double> lambda2, ArrayList<Integer> count2){
		
		int n=count2.size()-1;
		for (int i=0;i<3;i++){
			boolean found=false;
			for (;n>=0;n--){
				if(count2.get(n)>0){
					count2.set(n, 0);
					found=true;
					break;
				}
			}
			if(!found){
				break;
			}
		}
		
		return sign(lambda2,count2);
	}*/
	
	
	static Comparator<double[]> comp=(double[] a, double[] b)->{
		return -new Double(a[0]).compareTo(new Double(b[0]));
	};
	
	public static ArrayList<double[]> sort (ArrayList<Double> x, ArrayList<Integer> y){
		
		ArrayList<double[]> z=new ArrayList<double[]>();
		for (int i=0;i<x.size();i++){
			z.add(new double[]{x.get(i),y.get(i)});
		}
		Collections.sort(z,comp);
		return z;
	}
	//TODO: determine number of error mutations
	
	public static int sum(ArrayList<Integer> a){
		int sum=0;
		for (int i=0;i<a.size();i++){
			sum+=a.get(i);
		}
		return sum;
	}
	
	public static double sum_double(ArrayList<Double> a){
		double sum=0;
		for (int i=0;i<a.size();i++){
			sum+=a.get(i);
		}
		return sum;
	}
	
	
	
	
	
	
	public static double p(ArrayList<Double> lambda2, ArrayList<Integer> count2){
		
		
		int NN=0;
		double sum_lambda=0;
		for (int i=0;i<count2.size();i++){
			NN+=count2.get(i);
			sum_lambda+=lambda2.get(i);
		}
		if(lambda2.size()==0||NN==0){
			return 1;
		}
		
		double T=0;
		for (int i=0;i<count2.size();i++){
			if(count2.get(i)>0){
				T+=count2.get(i)*Math.log(lambda2.get(i));
				//System.out.println(Math.log(lambda2.get(i)));
			}
			if(count2.get(i)>0){
				T-=Gamma.logGamma(count2.get(i)+1);
			}
		}
		double offset=-NN*Math.log(sum_lambda);
		T+=offset;
		return T;
	}
	
	
	
	/*
	public static double sign_long(ArrayList<Double> lambda2, ArrayList<Integer> count2, double rate, double cov){
		
		
		int NN=0;
		for (int i=0;i<count2.size();i++){
			NN+=count2.get(i);
		}
		if(lambda2.size()==0||NN==0){
			return 1;
		}
		
		double T=0;
		for (int i=0;i<count2.size();i++){
			if(count2.get(i)>0){
				T+=count2.get(i)*Math.log(lambda2.get(i));
				//System.out.println(Math.log(lambda2.get(i)));
			}
			if(count2.get(i)>0){
				T-=Gamma.logGamma(count2.get(i)+1);
			}
		}
		double offset=NN*Math.log((double)(1)/(double)(lambda2.size()));
		T+=offset;
		//System.out.println(T);
		
		//System.out.println("T: "+Math.exp(T));
		Comparator<Double> comp=(Double a, Double b)->{
			return -a.compareTo(b);
		};
		//int NN=500;
		
		Collections.sort(lambda2,comp);
		double sum_lambda=0;
		for (int i=0;i<lambda2.size();i++){
			sum_lambda+=lambda2.get(i);
		}
		
		double sum_pow=0;
		{
			int k=2;
			for (int i=0;i<lambda2.size();i++){
				sum_pow+=Math.pow(lambda2.get(i),k);
			}
		}
	
		

		double fraction_high=0;
		ArrayList<Double> lambda_high=new ArrayList<Double>();
		{
			int k=2;
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),k);
				if(sum_local/sum_pow>0.99){
					
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
			}
		}
	
		double[] lambda_dist=prob_lambda(lambda2);
		
		
		double[] cum_high=new double[lambda_high.size()];
		double pp=0;
		for (int i=0;i<cum_high.length;i++){
			cum_high[i]=pp;
			pp+=lambda_high.get(i);
			
		}
		for (int i=0;i<cum_high.length;i++){
			cum_high[i]/=pp;
		}
		
		double[] cum_all=new double[lambda2.size()];
		pp=0;
		for (int i=0;i<cum_all.length;i++){
			cum_all[i]=pp;
			pp+=lambda2.get(i);
			
		}
		for (int i=0;i<cum_all.length;i++){
			cum_all[i]/=pp;
		}
	
		ArrayList<Double> s1=new ArrayList<Double>();
		ArrayList<Double> s3=new ArrayList<Double>();
		
		
		BinomialDistribution binom=new BinomialDistribution(NN,fraction_high);
		//System.out.println("SART");
		for (int k1=0;k1<10000;k1++){//10000
			int N_high=binom.sample();
			int[] meta_random=random_collision(N_high,cum_high);
			meta_random[0]+=NN-N_high;
			double sum1=0;
			for (int i=1;i<meta_random.length;i++){
				if(meta_random[i]>0){
					sum1+=Gamma.logGamma(i+1+1)*meta_random[i];
				}
			}
			s1.add(-sum1);
			//for (int i=0;i<meta_random.length;i++){
			//	System.out.print("	"+meta_random[i]);
			//}
			//System.out.println();
		}
		//System.out.println("END");
		
	
		
		PoissonDistribution dist=new PoissonDistribution (rate*cov);
		for (int k2=0;k2<100000;k2++){
			int err=Math.min(NN, dist.sample());
			s3.add(offset+random_draw(NN-err,lambda_dist,err,lambda2));
		}
		
		
		//System.out.println("END2");
		Collections.sort(s1);
		Collections.sort(s3);
		//System.out.println("END3");
		//System.out.println("Compare");
		
		
		
			int nn=0;
			for (int i=0;i<s3.size();i++){
				if(s3.get(i)+s1.get(0)>T){
					break;
				}
				for (int j=0;j<s1.size();j++){
					if(s3.get(i)+s1.get(j)>T){
						break;
					}
					nn++;
				}
			}
			return (double)(nn)/(double)(s1.size()*s3.size());
		

		
	}
	*/
	
	/*
	public static double sign_quick(ArrayList<Double> lambda2, ArrayList<Integer> count2, double rate, double cov){
		
		
		int NN=0;
		for (int i=0;i<count2.size();i++){
			NN+=count2.get(i);
		}
		if(lambda2.size()==0||NN==0){
			return 1;
		}
		
		double T=0;
		for (int i=0;i<count2.size();i++){
			if(count2.get(i)>0){
				T+=count2.get(i)*Math.log(lambda2.get(i));
				//System.out.println(Math.log(lambda2.get(i)));
			}
			if(count2.get(i)>0){
				T-=Gamma.logGamma(count2.get(i)+1);
			}
		}
		double offset=NN*Math.log((double)(1)/(double)(lambda2.size()));
		T+=offset;
		
		Comparator<Double> comp=(Double a, Double b)->{
			return -a.compareTo(b);
		};
		//int NN=500;
		
		Collections.sort(lambda2,comp);
		double sum_lambda=0;
		for (int i=0;i<lambda2.size();i++){
			sum_lambda+=lambda2.get(i);
		}
		
		double sum_pow=0;
		{
			int k=2;
			for (int i=0;i<lambda2.size();i++){
				sum_pow+=Math.pow(lambda2.get(i),k);
			}
		}
	
		

		double fraction_high=0;
		ArrayList<Double> lambda_high=new ArrayList<Double>();
		{
			int k=2;
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),k);
				if(sum_local/sum_pow>0.99){
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
			}
		}
	
		double[] lambda_dist=prob_lambda(lambda2);
		
		
		double[] cum_high=new double[lambda_high.size()];
		double pp=0;
		for (int i=0;i<cum_high.length;i++){
			cum_high[i]=pp;
			pp+=lambda_high.get(i);
			
		}
		for (int i=0;i<cum_high.length;i++){
			cum_high[i]/=pp;
		}
		
		double[] cum_all=new double[lambda2.size()];
		pp=0;
		for (int i=0;i<cum_all.length;i++){
			cum_all[i]=pp;
			pp+=lambda2.get(i);
			
		}
		for (int i=0;i<cum_all.length;i++){
			cum_all[i]/=pp;
		}
	
		ArrayList<Double> s1=new ArrayList<Double>();
		ArrayList<Double> s3=new ArrayList<Double>();
		
		
		BinomialDistribution binom=new BinomialDistribution(NN,fraction_high);
		//System.out.println("SART");
		for (int k1=0;k1<100;k1++){//
			int N_high=binom.sample();
			int[] meta_random=random_collision(N_high,cum_high);
			meta_random[0]+=NN-N_high;
			double sum1=0;
			for (int i=1;i<meta_random.length;i++){
				if(meta_random[i]>0){
					sum1+=Gamma.logGamma(i+1+1)*meta_random[i];
				}
			}
			s1.add(-sum1);
			//for (int i=0;i<meta_random.length;i++){
			//	System.out.print("	"+meta_random[i]);
			//}
			//System.out.println();
		}
		//System.out.println("END");
		
		PoissonDistribution dist=new PoissonDistribution (rate*cov);
		for (int k2=0;k2<1000;k2++){
			int err=Math.min(NN, dist.sample());
			s3.add(offset+random_draw(NN-err,lambda_dist,err,lambda2));
		}
		
		Collections.sort(s1);
		Collections.sort(s3);
		
		int n1=0;
		int n2=0;
		for (int k=0;k<1000;k++){
			if(s3.get((int)(s3.size()*Math.rrandom()))+s1.get((int)(s1.size()*Math.rrandom()))<=T){
				n1++;
			}
			n2++;
		}
		
		return (double)(n1)/(double)(n2);
			
		
	}
	*/
	
	public static  int overlap (String[] tt, ArrayList<Integer> t){
		int n=0;
		int[] ttt=new int[tt.length];
		for (int i=0;i<tt.length;i++){
			ttt[i]=Integer.parseInt(tt[i]);
		}
		
		for (int i=0;i<t.size();i++){
			for (int j=0;j<ttt.length;j++){
				if(t.get(i).intValue()==ttt[j]){
					n++;
				}
			}
		
		}	
		
		return n;
	}
	
	
	
	
	
	
	
	public static double[] prob_lambda(ArrayList<Double> lambda_high){
		if(lambda_high.size()==0){
			return new double[0];
		}
		double sum=0;
		for (int i=0;i<lambda_high.size();i++){
			sum+=lambda_high.get(i);
		}
		
		double[] array=new double[100000];
		for (int i=0;i<array.length;i++){
			array[i]=Double.NaN;
		}
		double cum=0;
		for (int i=0;i<lambda_high.size();i++){
			array[(int)(cum*(array.length-1))]=Math.log(lambda_high.get(i));
			cum+=lambda_high.get(i)/sum;
		}
		
		double prev=Double.NaN;
		for (int i=array.length-1;i>=0;i--){
			if(Double.isNaN(array[i])){
				array[i]=prev;
			}
			else{
				prev=array[i];
			}
		}
		int i=array.length-1;
		while(i>=0&&Double.isNaN(array[i])){
			i--;
		}
		for (int ii=i+1;ii<array.length;ii++){
			array[ii]=array[i];
		}
		
		return array;
		
		
	}
	
	/*
	public static double random_draw(int N, double[] lambda){
		if(N==0){
			return 0;
		}
		double sum=0;
		for (int i=0;i<N;i++){
			sum+=lambda[(int)(lambda.length*Math.rrandom())];
		}
		return sum;
	}*/
	
	
	
	/*
	public static double random_draw(int N, double[] lambda, double rate, double cov, ArrayList<Double> lambda2){
		if(N==0){
			return 0;
		}
		double sum=0;
		for (int i=0;i<N;i++){
			sum+=lambda[(int)(lambda.length*Math.random())];
		}
		
		
		for (int i=0;i<N_err;i++){
			sum+=Math.log(lambda2.get((int)(Math.random()*lambda2.size())));
		}
		return sum;
	}
	*/
	
	public static int[] random_meta(BinomialDistribution[] p_meta_cum, int N_high){
		int[] count=new int[p_meta_cum.length];
		for (int i=p_meta_cum.length-1;i>0;i--){
			int a=p_meta_cum[i].sample();
			while(a*(i+1)>N_high){
				System.out.println("AR");
				 a=p_meta_cum[i].sample();
			}
			count[i]=a;
			N_high-=a*(i+1);
		}
		count[0]=N_high;
		return count;
	}
	
	/*public static int[] random_meta(double[] p_meta_cum, int N_high){
		int[] count_power=new int[p_meta_cum.length];
		 for (int i=0;i<N_high;){
			 double r=Math.random();
			 for (int j=0;j<p_meta_cum.length;j++){
				 if(j<p_meta_cum.length-1){
						if(p_meta_cum[j]<=r&&r<p_meta_cum[j+1]){
							if(i+(j+1)<=N_high){
								count_power[j]++;
								i+=j+1;
							}
							break;
						}
				}
					else{
						if(p_meta_cum[p_meta_cum.length-1]<=r){
							if(i+(j+1)<N_high){
								count_power[count_power.length-1]++;
								i+=j+1;
							}
							break;
						}
					}	
			}	
		
		 }
		return count_power;	
	}*/
	
	public static double[] prob(double[] power, int N){
		double[] p=new double[power.length];
		for (int k=1;k<=power.length;k++){
			p[k-1]=power[k-1]*Math.exp(Gamma.logGamma(N+1)-Gamma.logGamma(k+1)-Gamma.logGamma(N-k+1));
		}
		double sum=0;
		for (int i=0;i<p.length;i++){
			sum+=p[i];
		}
		for (int i=0;i<p.length;i++){
			p[i]/=sum;
		}
		return p;
		
	}
	
	public static double[] cum(double[] p){
		double cum=0;
		double[] pp=new double[p.length];
		for (int i=0;i<p.length;i++){
			pp[i]=cum;
			cum+=p[i];
		}
		for (int i=0;i<p.length;i++){
			pp[i]/=cum;
		}
		return pp;
	}
	

	public static int[] meta(int[] v){
		int[] meta=new int[6];
		for (int i=0;i<v.length;i++){
			if(v[i]<meta.length){
				meta[v[i]]++;
			}
		}
		return meta;
	}
	
	public static double[] prob(ArrayList<Double> prob){
		double sum=0;
		for (int i=0;i<prob.size();i++){
			sum+=prob.get(i);
		}
		ArrayList<Double> prob2=new ArrayList<Double>();
		for (int i=0;i<prob.size();i++){
			prob2.add(prob.get(i)/sum);
		}
		Collections.sort(prob2);
		double[] array=new double[100000];
		for (int i=0;i<array.length;i++){
			array[i]=Double.NaN;
		}
		double cum=0;
		for (int i=0;i<prob2.size();i++){
			//if(Double.isNaN(array[(int)(cum*prob2.size())])){
				array[(int)(cum*array.length)]=Math.log(prob2.get(i));
			//}
			cum+=prob2.get(i);
		}
		
		double prev=Double.NaN;
		for (int i=array.length-1;i>=0;i--){
			if(Double.isNaN(array[i])){
				array[i]=prev;
			}
			else{
				prev=array[i];
			}
		}
		int i=array.length-1;
		while(Double.isNaN(array[i])){
			i--;
		}
		for (int ii=i+1;ii<array.length;ii++){
			array[ii]=array[i];
		}
		
		return array;
		
	}
}

