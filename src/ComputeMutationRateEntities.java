/************************************************************           
 * MutPanning - Step 7										*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2018 									*   
 *															*   
 * License:		Public Domain								*   
 *															*   
 * Summary: This script fulfills two 2 purposes:			*
 * A) It computes the local mutation rate of each cancer 	*
 * entity. For this purpose, it first computes the fraction	*
 * of mutations that are contributed to the no. mutations	*
 * which are contributed from each cluster to the total no.	*
 * mutations for each cancer entities. Then the script sums	*
 * the mutation rates up according to these weights which	*
 * delivers the total mutation rate of each cluster.		*
 * B) To calibrate the CBASE model which determines the 	*
 * statistics whether a mutation contains more mutations	*
 * than expected, we needed for each gene a "target size"	*
 * as well as its mutation count, separately for syn. and	*
 * nonsynonymous mutations. While the script walks through	*
 * the background mutation rates, it computes these values	*
 * for each gene on the go.									*
 *  Note that this script provides the files that will be	*
 * immediately for the statistics of the mutation counts	*
 * (CBASE model) as well as the sequence context			*
 * 															*   
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


public class ComputeMutationRateEntities {
	static String[] chr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};	
	
	
	static String file_samples="";
	static String file_align="";
	static String file_genes="";
	static String file_type="";
	static String file_clusters="";
	static String file_lambda="";
	static String file_out="";
	static String file_out2="";
	
	static String[] entities=new String[0];
	static String[] index_header_samples={"ID","Sample","Cohort"};
	
	
	static double[][] weights=new double[0][0];
	static double coverage_threshold=0.05;
	static int no_clusters=-1;
	
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
			"ASMT",//double genes
			"ASMTL",
			"CSF2RA",
			"DHRSX",
			"IL3RA",
			"IL3R",
			"IL9R",
			"SPRY3",
			"ZBED1"
		};
	
	
	static ArrayList<Gene>[] genes=new ArrayList[chr.length];
	static ArrayList<Integer>[] list=new ArrayList[0];
	
	/*
	 * argument0: root file
	 * argugment1: sample file
	* argument2: path to the Hg19 folder
	 */
	
	public static void main(String[] args){
		
		file_samples=args[1];
		file_align=args[0]+"AlignHg19/AlignHg19Chr";
		file_genes=args[2]+"Exons_Hg19.txt";
		file_type=args[0]+"AffinityCounts/TypeCount.txt";
		file_clusters=args[0]+"ClusteringComplete/ClusteringComplete_Samples.txt";
		file_lambda=args[0]+"MutationRateClusters/Lambda_Chr";
		file_out=args[0]+"EntityCounts/EntityCounts";
		file_out2=args[0]+"CBASE/CountsRaw/Count";
		
		if (!new File(args[0]+"EntityCounts/").exists()){
			new File(args[0]+"EntityCounts/").mkdir();
		}
		if (!new File(args[0]+"CBASE/CountsRaw/").exists()){
			new File(args[0]+"CBASE/CountsRaw/").mkdirs();
		}
		
		//Determine all entity names, as clustering is performed separately for each Step this is needed to coordinate the order
		try{
			FileInputStream in=new FileInputStream(file_samples);
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
			entities=new String[aa.size()];
			for (int i=0;i<aa.size();i++){
				entities[i]=aa.get(i);
			}
		}
		catch(Exception e){
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			System.out.println(e);
		}		
		
		
		try{
			//this part is to determine how much which cluster relatively
			//contributes to the mutation rate of the entity. for this purpose
			//the no. mutations are summed up in each cluster and compared to
			//the total number of mutations in the entity
			String s="";
			ArrayList<String> names_count=new ArrayList<String>(); 
			ArrayList<Integer> count=new ArrayList<Integer>();
			
			//read the no. mutaions per samples
			FileInputStream in=new FileInputStream(file_type);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				names_count.add(t[0]);
				int c=0;
				for (int i=1;i<=6;i++){
					c+=Integer.parseInt(t[i]);
				}
				count.add(c);
			}
			input.close();
			
			no_clusters=0;
			in=new FileInputStream(file_clusters);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			while((s=input.readLine())!=null){
				if(Integer.parseInt(s.split("	")[2])>no_clusters){
					no_clusters=Integer.parseInt(s.split("	")[2]);
				}
			}
			no_clusters++;
			input.close();
			
			weights=new double[entities.length+1][no_clusters];
			
			//sum up the no. mutations in each cluster
			System.out.println("START");
			in=new FileInputStream(file_clusters);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				int c=count.get(index(t[0],names_count));
				int index_entity=index(t[1],entities);
				int cluster=Integer.parseInt(t[2]);
				
				weights[index_entity][cluster]+=c;
				weights[weights.length-1][cluster]+=c;
			}
			input.close();
			
			
			//normalize weights. weight <0.01 are set to 0
			//this safes a lot of run time and does not affect the
			//mutation rate of the sample substantially
			for (int i=0;i<weights.length;i++){
				double sum=0;
				for (int j=0;j<weights[i].length;j++){
					sum+=weights[i][j];
				}
				for (int j=0;j<weights[i].length;j++){
					weights[i][j]/=sum;
				}
				for (int j=0;j<weights[i].length;j++){
					if(weights[i][j]<0.01){
						weights[i][j]=0;
					}
				}
				sum=0;
				for (int j=0;j<weights[i].length;j++){
					sum+=weights[i][j];
				}
				for (int j=0;j<weights[i].length;j++){
					weights[i][j]/=sum;
				}
			}
			
			
			for (int i=0;i<genes.length;i++){
				genes[i]=new ArrayList<Gene>();
			}
				
			
			//read genes with positions
			in=new FileInputStream(file_genes);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			
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
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					genes[i].get(j).start=min(genes[i].get(j).coord);
					genes[i].get(j).end=max(genes[i].get(j).coord);
				}
			}
			
			list=new ArrayList[entities.length];
			for (int i=0;i<list.length;i++){
				list[i]=new ArrayList<Integer>();
			}
			in=new FileInputStream(file_samples);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				int ii=index(t[index_header[2]],entities);
				list[ii].add(Integer.parseInt(t[index_header[0]]));
			}
			input.close();
			
			
			//run gene count / mutation rate computation on each chr independently
			Subthread[] threads=new Subthread[chr.length];
			for (int i=0;i<threads.length;i++){
				threads[i]=new Subthread();
				threads[i].c=i;
				threads[i].start();
			}

			//wait until all threads are finished before closing the scripts
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
			
			
		}
		catch(Exception e){
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			System.out.println(e);
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
	
	public static boolean contains(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
	
	public static int index(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
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
	public static int index_gene(String s, ArrayList<Gene> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).name.equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	
	//each thread handles the counts of the genes & computation of the 
	//entity specific mutation rate independently
	private static class Subthread extends Thread{
		volatile boolean done=false;
		int c=-1;
		
		public void run(){
			System.out.println("start "+chr[c]);
			try{
				ArrayList<Integer> position=new ArrayList<Integer>();
				ArrayList<String> nucl=new ArrayList<String>();
				ArrayList<Double> coverage=new ArrayList<Double>();
				ArrayList<int[]> label=new ArrayList<int[]>();
				ArrayList<int[][]> count=new ArrayList<int[][]>();
				ArrayList<double[][]> lambda=new ArrayList<double[][]>();
				
				FileInputStream in=new FileInputStream(file_lambda+chr[c]+".txt");
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				
				
				FileInputStream in2=new FileInputStream(file_align+chr[c]+".txt");
				DataInputStream inn2=new DataInputStream(in2);
				BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
				
				FileWriter[] out=new FileWriter[entities.length+1];
				BufferedWriter[] output= new BufferedWriter[entities.length+1];
				
				FileWriter[] out2=new FileWriter[entities.length+1];
				BufferedWriter[] output2= new BufferedWriter[entities.length+1];
				
				for (int i=0;i<entities.length;i++){
					out[i]=new FileWriter(file_out+entities[i]+"_Chr"+chr[c]+".txt");
					output[i]=new BufferedWriter(out[i]);
					out2[i]=new FileWriter(file_out2+entities[i]+"_Chr"+chr[c]+".txt");
					output2[i]=new BufferedWriter(out2[i]);
				}
				out[out.length-1]=new FileWriter(file_out+"PanCancer"+"_Chr"+chr[c]+".txt"); 
				output[out.length-1]=new BufferedWriter(out[out.length-1]);
				out2[out.length-1]=new FileWriter(file_out2+"PanCancer"+"_Chr"+chr[c]+".txt"); 
				output2[out.length-1]=new BufferedWriter(out2[out2.length-1]);
				
				
				//walk through the reference sequence of each gene
				for (int nn=0;nn<genes[c].size();nn++){
					int ii=0;
					if(position.size()>0){
						while(ii<position.size()&&position.get(ii)<genes[c].get(nn).start){
							ii++;
						}
					}
					
					//remove the squence of the previous gene out of the queue
					if(ii>0){
						for (int i=ii-1;i>=0;i--){
							nucl.remove(i);
							lambda.remove(i);
							position.remove(i);
							coverage.remove(i);
							count.remove(i);
							label.remove(i);
						}
					}
					
					//read the gene reference sequence together with the clusterwise mutation rate
					String s2="";
					while ((s2=input2.readLine())!=null){
						String s=input.readLine();
						//System.out.println(s);
						String[] t1=s.split("	");
						double[][] ll=new double[no_clusters][3];
						for (int i=0;i<ll.length;i++){
							String[] tt=t1[i+2].split(";");
							for (int j=0;j<tt.length;j++){
								ll[i][j]=Double.parseDouble(tt[j]);
							}
						}
						
						String[] t2=s2.split("	");
						if(!t1[0].equals(t2[0])||!t1[1].equals(t2[1])){
							System.out.println("S: "+s);
							System.out.println("S2: "+s);
							System.out.println("0: "+t1[0].equals(t2[0]));
							System.out.println("1: "+t1[1].equals(t2[1]));
							System.out.println(t1[0]);
							System.out.println(t2[0]);
							System.out.println(t1[1]);
							System.out.println(t2[1]);
							System.out.println("incompatible");
							System.exit(0);
						}
						
						int label1=-1;
						int label2=-1;
						int label3=-1;
						
						if(t2.length>6){
							if(t2[6].equals(t2[8])){
								label1=0;//"syn";
							}
							else{
								label1=1;//"ns";
							}
						}
						else{
							label1=2;//"nc";
						}
						
						if(t2.length>6){
							if(t2[6].equals(t2[9])){
								label2=0;//"syn";
							}
							else{
								label2=1;//"ns";
							}
						}
						else{
							label2=2;//"nc";
						}
						
						if(t2.length>6){
							if(t2[6].equals(t2[10])){
								label3=0;//"syn";
							}
							else{
								label3=1;//"ns";
							}
						}
						else{
							label3=2;//"nc";
						}
						
						int[] label_local=new int[]{label1,label2,label3};
						
						int[] index1=new int[0];
						if(t2.length>3&&!t2[3].equals("")){
							index1=integer(t2[3].split(";"));
						}
								
						int[] index2=new int[0];
						if(t2.length>4&&!t2[4].equals("")){
							index2=integer(t2[4].split(";"));
						}
						int[] index3=new int[0];
						if(t2.length>5&&!t2[5].equals("")){
							index3=integer(t2[5].split(";"));
						}
						int[][] indices=new int[][]{index1,index2,index3};
						
						
						
						//combine the mutation rates according to the weights for each cancer entity
						double[][] lambda_local=new double[entities.length+1][3];
						for (int k=0;k<3;k++){
							for (int i=0;i<weights.length;i++){//cohorts
								for (int j=0;j<weights[i].length;j++){//clusters
									lambda_local[i][k]+=weights[i][j]*ll[j][k];
								}
							}
						}
						
						//count the no mutant samples for each entity, the last index refers
						//to the pan cancer count (i.e. all samples)
						int[][] count_local=new int[entities.length+1][3];
						for (int k=0;k<indices.length;k++){
							for (int i=0;i<list.length;i++){
								count_local[i][k]=overlap(indices[k],list[i]);
							}
							count_local[count_local.length-1][k]=indices[k].length;
						}
						
						//safe local mutation rate in queue, this will be later needed
						//to compute the target sizes and mutation count of gene
						position.add(Integer.parseInt(t2[0]));
						nucl.add(t2[1]);
						coverage.add(Double.parseDouble(t2[2]));
						label.add(label_local);
						count.add(count_local);
						lambda.add(lambda_local);
						//output of the local mutation rate && local mutation count for each
						//cancer entity
						for (int k=0;k<output.length;k++){
							output[k].write(t2[0]+"	"+t2[1]+"	"+t2[2]+"	"+ou(label_local)+"	"+ou(lambda_local[k])+"	"+ou(count_local[k]));
							output[k].newLine();
						}
						
						if(Integer.parseInt(t2[0])>=genes[c].get(nn).end){
							break;
						}
					}
					
					double cov_gene=0;
					double cov_gene_syn=0;
					int[] count_gene=new int[entities.length+1];
					int[] count_gene_syn=new int[entities.length+1];
					double[] lambda_gene=new double[entities.length+1];
					double[] lambda_gene_syn=new double[entities.length+1];
					
					//sum up the mutation rates and mutation counts over the gene to
					//compute its target size and total no. mutations
					//note that the mutation rate and target size is determined separately
					//for synonymous and nonsynonymous mutations
					for (int i=0;i<position.size();i++){
						if(coverage.get(i)<coverage_threshold){
						//	continue;
						}
						if(genes[c].get(nn).contains(position.get(i))){
							for (int j=0;j<label.get(i).length;j++){
								if(label.get(i)[j]!=0){
									cov_gene_syn+=coverage.get(i);
									for (int k=0;k<count.get(i).length;k++){
										count_gene[k]+=count.get(i)[k][j];
										lambda_gene[k]+=coverage.get(i)*lambda.get(i)[k][j];
									}
								}
								else{
									cov_gene+=coverage.get(i);
									for (int k=0;k<count.get(i).length;k++){
										count_gene_syn[k]+=count.get(i)[k][j];
										lambda_gene_syn[k]+=coverage.get(i)*lambda.get(i)[k][j];
									}
								}
							}
						}
					}
					
					//output of the gene count and target size
					//genes which have nearly no coverage are masked as they are producing artifacts
					if(cov_gene+cov_gene_syn>30){
						for (int i=0;i<output2.length;i++){
							output2[i].write(genes[c].get(nn).name+"	"+cov_gene+"	"+cov_gene_syn+"	"+lambda_gene[i]+"	"+lambda_gene_syn[i]+"	"+count_gene[i]+"	"+count_gene_syn[i]);
							output2[i].newLine();
						}
					}
					
					
				}
				input.close();
				input2.close();
				for (int i=0;i<output.length;i++){
					output[i].close();
					output2[i].close();
				}
				
				
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
	}
	
	public static int overlap(int[] a, ArrayList<Integer> b){
		int n=0;
		for (int i=0;i<a.length;i++){
			for (int j=0;j<b.size();j++){
				if(a[i]==b.get(j)){
					n++;
				}
			}
		}
		return n;
	}
	
	public static String ou (double[] a){
		if(a.length==0){
			return "";
		}
		if(a.length==1){
			return ""+a[0];
		}
		String s=""+a[0];
		for (int i=1;i<a.length;i++){
			s=s+";"+a[i];
		}
		return s;
	}
	
	public static String ou(int[] a){
		if(a.length==0){
			return "";
		}
		if(a.length==1){
			return ""+a[0];
		}
		String s=""+a[0];
		for (int i=1;i<a.length;i++){
			s=s+";"+a[i];
		}
		return s;
	}
	
	public static int[] integer(String[] s){
		int[] a=new int[s.length];
		for (int i=0;i<a.length;i++){
			a[i]=Integer.parseInt(s[i]);
		}
		return a;
	}
	
	
	public static int index(String s, String[] t){
		for (int i=0;i<t.length;i++){
			if(t[i].equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	private static class Gene{
		String name="";
		ArrayList<int[]> coord=new ArrayList<int[]>();
		int start=-1;
		int end=-1;
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
}
