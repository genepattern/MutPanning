/************************************************************           
 * MutPanning - Step 6									*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2018 									*   
 *															*   
 * License:		Public Domain								*   
 *															*   
 * Summary: This script computes the mutation rate of each	*
 * 			cluster. These clusters were dervied in the 	*
 * 			previous steps and contain samples with 		*
 * 			similar passenger mutation distribution patterns*
 * 			Mutation rates are computed for all clusters in	*
 * 			parallel as a major timelimiting factor is 		*
 *			reading the sequence context around each		*
 *			position. Each chr is handled in a separate 	*
 *			thread.											*
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


public class ComputeMutationRateClusters {
	static ArrayList<double[][][]> lambda_context=new ArrayList<double[][][]>();
	static ArrayList<double[][]> lambda_type=new ArrayList<double[][]>();
	static String[] chr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
	
	
	static String file_out="";
	static String file_seq="";
	static String file_signatures="";
	static String file_reference="";
	
	/*
	 * argument0: root file
	 * argument1: path to the Hg19 folder
	 */
	
	public static void main(String[] args){
		
		file_out=args[0]+"MutationRateClusters/Lambda_Chr";
		file_seq=args[1]+"ASAnnotationHg19/ASAnnotation_chr";
		file_signatures=args[0]+"ClusteringComplete/ClusteringComplete_Affinity.txt";
		file_reference=args[1]+"FileReferenceCount.txt";
		
		
		try{
			if(!new File(args[0]+"MutationRateClusters/").exists()){
				new File(args[0]+"MutationRateClusters/").mkdir();
			}
			
			//each cluster delivers a "signature" for the distribution of passenger
			//mutations in that cluster. read that signature
			ArrayList<double[][][]> frequency=new ArrayList<double[][][]>();
			FileInputStream in=new FileInputStream(file_signatures);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			String s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				int[][][] a =new int [20][6][4];
				int n=0;
				for (int i=0;i<a.length;i++){
					for (int j=0;j<a[i].length;j++){
						for (int k=0;k<a[i][j].length;k++){
							a[i][j][k]=Integer.parseInt(t[n+1]);
							n++;
						}
					}
				}
				int[] b=new int[6];
				int sum=0;
				for (int i=0;i<6;i++){
					for (int j=0;j<4;j++){
						b[i]+=a[9][i][j]+a[10][i][j];
						sum+=a[9][i][j]+a[10][i][j];
						//System.out.println(sum);
					}
				}
				frequency.add(freq(a));
				double[][] lambda_t=new double[2][];
				lambda_t[0]=new double[]{2*(double)(b[0]+b[1]+b[2])/(double)(sum),2*(double)(b[3]+b[4]+b[5])/(double)(sum)};
				lambda_t[1]=new double[]{3*(double)(b[0]+b[3])/(double)(sum),3*(double)(b[1]+b[5])/(double)(sum),3*(double)(b[2]+b[4])/(double)(sum)};
				lambda_type.add(lambda_t);
			}
			input.close();
			
			
			//read the distribution in the reference genome
			int[][][] reference_abs=new int[20][2][4];
			in=new FileInputStream(file_reference);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			for (int i=0;i<2;i++){
				for (int j=0;j<4;j++){
					String[] t=input.readLine().split("	");
					for (int k=0;k<20;k++){
						if(k<10){
							reference_abs[k][i][j]=Integer.parseInt(t[k+2]);
						}
						else{
							reference_abs[k][i][j]=Integer.parseInt(t[k+3]);
						}
					}
				}
			}
			input.close();
			
			//Compute the likelihood ratios lambda as the ratio between reference and observed counts
			double[][][] reference=freq(reference_abs);
			for (int i=0;i<frequency.size();i++){
				double[][][] ll=new double[20][6][4];
				for (int j=0;j<frequency.get(i).length;j++){
					for (int k=0;k<frequency.get(i)[j].length;k++){
						for (int l=0;l<frequency.get(i)[j][k].length;l++){
							ll[j][k][l]=frequency.get(i)[j][k][l]/reference[j][k/3][l];
						}
					}
				}
				lambda_context.add(ll);
			}
			
			
			
			//compute the distribution at each position for each chr separately
			Subthread[] threads=new Subthread[chr.length];
			for (int i=0;i<threads.length;i++){
				threads[i]=new Subthread();
				threads[i].c=i;
				threads[i].start();
			}
			
			//wait until done for all chr before closing the script
			boolean all_done=false;
			do{
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
	
	//each thread computes the clusterwise local mutation rate for 1 chromosome
	//(reading the reference sequence is usually the rate limiting step)
	private static class Subthread extends Thread{
		
		ArrayList<Integer> pos=new ArrayList<Integer>();
		ArrayList<String> nucl=new ArrayList<String>();
		int c=-1;
		BufferedWriter output=null;
		volatile boolean done=false;
		public void run(){
			try{
				String s="";
				for (int i=0;i<50;i++){
					pos.add(-1);
					nucl.add("N");
				}
				
				//walk through the refernce sequence and safe it into a queue
				FileWriter out=new FileWriter(file_out+chr[c]+".txt");
				output= new BufferedWriter(out);
				FileInputStream in=new FileInputStream(file_seq+chr[c]+".txt");
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					pos.add(Integer.parseInt(t[0]));
					nucl.add(t[1]);
					
					//when the queue is large enough, compute the mutation rate in the cnenter
					//(update function) and delete the first element
					if(pos.size()>100){
						update();
						pos.remove(0);
						nucl.remove(0);
					}
					
				}
				
				//compute the last mutation rates at the end of the chr until the queue empty 
				for (int i=0;i<50;i++){
					pos.add(-1);
					nucl.add("N");
					update();
					pos.remove(0);
					nucl.remove(0);
				}
				input.close();
				output.close();
			}
			catch(Exception e){
				StackTraceElement[] aa=e.getStackTrace();
				for (int i=0;i<aa.length;i++){
					System.out.println(i+"	"+aa[i].getLineNumber());
				}
				System.out.println(e);
			}
			done=true;
			System.out.println("Done "+chr[c]);
		}
		
		
		//the "heart" of this thread which computes the local mutation rate in the center of a queu and appends it to output
		//in brief walks from -10 to +10 around the central position (skipping 0) and forming the product over likelihoods
		public void update(){
			try{
				if(nucl.get(50).equals("A")){
					double[][] lambda=new double[lambda_context.size()][3];
					for (int i=0;i<lambda.length;i++){
						lambda[i][0]=lambda_type.get(i)[1][0];
						lambda[i][1]=lambda_type.get(i)[1][1];
						lambda[i][2]=lambda_type.get(i)[1][2];
						
						lambda[i][0]*=lambda_type.get(i)[0][1];
						lambda[i][1]*=lambda_type.get(i)[0][1];
						lambda[i][2]*=lambda_type.get(i)[0][1];
						int[] tt={3,5,4};
						for (int j=-10;j<0;j++){
							for (int k=0;k<tt.length;k++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									lambda[i][k]*=lambda_context.get(i)[9-j][tt[k]][3-nucl_index(nucl.get(50+j))];
								}
							}
						}
						for (int j=1;j<=10;j++){
							for (int k=0;k<tt.length;k++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									lambda[i][k]*=lambda_context.get(i)[10-j][tt[k]][3-nucl_index(nucl.get(50+j))];
								}
							}
						}
					}
					for (int k=0;k<lambda.length;k++){
						for (int l=0;l<lambda[k].length;l++){
							if(lambda[k][l]>1000){
								lambda[k][l]=1000;
							}
						}
					}
					output.write(pos.get(50)+"	"+nucl.get(50));
					for (int k=0;k<lambda.length;k++){
						output.write("	"+lambda[k][0]+";"+lambda[k][1]+";"+lambda[k][2]);
					}
					output.newLine();
				}
				else if(nucl.get(50).equals("C")){
					double[][] lambda=new double[lambda_context.size()][3];
					for (int i=0;i<lambda.length;i++){
						lambda[i][0]=lambda_type.get(i)[1][0];
						lambda[i][1]=lambda_type.get(i)[1][1];
						lambda[i][2]=lambda_type.get(i)[1][2];
						
						lambda[i][0]*=lambda_type.get(i)[0][0];
						lambda[i][1]*=lambda_type.get(i)[0][0];
						lambda[i][2]*=lambda_type.get(i)[0][0];
						int[] tt={0,1,2};
						for (int j=-10;j<0;j++){
							for (int k=0;k<tt.length;k++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									lambda[i][k]*=lambda_context.get(i)[10+j][tt[k]][nucl_index(nucl.get(50+j))];
								}
							}
						}
						for (int j=1;j<=10;j++){
							for (int k=0;k<tt.length;k++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									lambda[i][k]*=lambda_context.get(i)[9+j][tt[k]][nucl_index(nucl.get(50+j))];
								}
							}
						}
					}
					for (int k=0;k<lambda.length;k++){
						for (int l=0;l<lambda[k].length;l++){
							if(lambda[k][l]>1000){
								lambda[k][l]=1000;
							}
						}
					}
					output.write(pos.get(50)+"	"+nucl.get(50));
					for (int k=0;k<lambda.length;k++){
						output.write("	"+lambda[k][0]+";"+lambda[k][1]+";"+lambda[k][2]);
					}
					output.newLine();
				}
				else if(nucl.get(50).equals("G")){
					double[][] lambda=new double[lambda_context.size()][3];
					for (int i=0;i<lambda.length;i++){
						lambda[i][0]=lambda_type.get(i)[1][0];
						lambda[i][1]=lambda_type.get(i)[1][1];
						lambda[i][2]=lambda_type.get(i)[1][2];
						
						lambda[i][0]*=lambda_type.get(i)[0][0];
						lambda[i][1]*=lambda_type.get(i)[0][0];
						lambda[i][2]*=lambda_type.get(i)[0][0];
						int[] tt={0,1,2};
						for (int j=-10;j<0;j++){
							for (int k=0;k<tt.length;k++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									lambda[i][k]*=lambda_context.get(i)[9-j][tt[k]][3-nucl_index(nucl.get(50+j))];
								}
							}
						}
						for (int j=1;j<=10;j++){
							for (int k=0;k<tt.length;k++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									lambda[i][k]*=lambda_context.get(i)[10-j][tt[k]][3-nucl_index(nucl.get(50+j))];
								}
							}
						}
					}
					for (int k=0;k<lambda.length;k++){
						for (int l=0;l<lambda[k].length;l++){
							if(lambda[k][l]>1000){
								lambda[k][l]=1000;
							}
						}
					}
					output.write(pos.get(50)+"	"+nucl.get(50));
					for (int k=0;k<lambda.length;k++){
						output.write("	"+lambda[k][0]+";"+lambda[k][1]+";"+lambda[k][2]);
					}
					output.newLine();
				}
				else if(nucl.get(50).equals("T")){
					double[][] lambda=new double[lambda_context.size()][3];
					for (int i=0;i<lambda.length;i++){
						lambda[i][0]=lambda_type.get(i)[1][0];
						lambda[i][1]=lambda_type.get(i)[1][1];
						lambda[i][2]=lambda_type.get(i)[1][2];
						
						lambda[i][0]*=lambda_type.get(i)[0][1];
						lambda[i][1]*=lambda_type.get(i)[0][1];
						lambda[i][2]*=lambda_type.get(i)[0][1];
						int[] tt={3,5,4};
						for (int j=-10;j<0;j++){
							for (int k=0;k<tt.length;k++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									lambda[i][k]*=lambda_context.get(i)[10+j][tt[k]][nucl_index(nucl.get(50+j))];
								}
							}
						}
						for (int j=1;j<=10;j++){
							for (int k=0;k<tt.length;k++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									lambda[i][k]*=lambda_context.get(i)[9+j][tt[k]][nucl_index(nucl.get(50+j))];
								}
							}
						}
					}
					for (int k=0;k<lambda.length;k++){
						for (int l=0;l<lambda[k].length;l++){
							if(lambda[k][l]>1000){
								lambda[k][l]=1000;
							}
						}
					}
					output.write(pos.get(50)+"	"+nucl.get(50));
					for (int k=0;k<lambda.length;k++){
						output.write("	"+lambda[k][0]+";"+lambda[k][1]+";"+lambda[k][2]);
					}
					output.newLine();
				}
				
			}
			catch(Exception e){
				StackTraceElement[] aa=e.getStackTrace();
				for (int i=0;i<aa.length;i++){
					System.out.println(i+"	"+aa[i].getLineNumber());
				}
				System.out.println(e);
			}
			
			
		}
		
		public int nucl_index(String s){
			if(s.toUpperCase().equals("A")){
				return 0;
			}
			else if(s.toUpperCase().equals("C")){
				return 1;
			}
			else if(s.toUpperCase().equals("G")){
				return 2;
			}
			else if(s.toUpperCase().equals("T")){
				return 3;
			}
			return -1;
		}
	}
	
	public static double[][][] freq(int[][][] a){
		double[][][] b=new double[a.length][a[0].length][a[0][0].length];
		for (int i=0;i<a.length;i++){
			for (int j=0;j<a[i].length;j++){
				int sum=0;
				for (int k=0;k<a[i][j].length;k++){
					sum+=a[i][j][k];
				}
				if(sum>=50){
					for (int k=0;k<a[i][j].length;k++){
						if(a[i][j][k]==0){
							a[i][j][k]++;
							sum++;
						}
					}
					for (int k=0;k<a[i][j].length;k++){
						b[i][j][k]=(double)(a[i][j][k])/(double)(sum);
					}
				}
				else{
					
					for (int k=0;k<a[i][j].length;k++){
						b[i][j][k]=(double)((double)(50-sum)/(double)(a[i][j].length)+a[i][j][k])/(double)(50);
								
					}
					
					
					/*for (int k=0;k<a[i][j].length;k++){
						b[i][j][k]=1.0/(double)(a[i][j].length);//(double)(a[i][j][k])/(double)(sum);
					}*/
				}
			}
		}
		return b;
	}
}
