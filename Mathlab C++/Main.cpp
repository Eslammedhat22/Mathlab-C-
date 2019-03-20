#include "CMatrix.h"
#include <iostream>
#include <string>
#include <cstring>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <cstdio>
#include <exception>
#include <stack>
#include <string.h>
#include "math.h"
#include <ctgmath>
using namespace std;

//global variables 
int NVar = 0;
int tempNVar = 0;
//array of varibles names
string * varibleNames = new string[100];
//array of Matrices
CMatrix * Matrices = new CMatrix[100];
const string  tempvaribleNames [20]= {"#","##","###","####","#####","######","#######","########","#########","##########","###########","############","#############","##############","##############",
                                      "################","#################","##################","###################","####################"};
//array of  temp Matrices
CMatrix * tempMatrices = new CMatrix[20];
const string functionsname []={ "sinh",	"cosh",	"tanh","coth","sech","csch",
	                            "sin","cos","tan","cot","sec","csc",	
                                "asin",	"acos","atan","acot","asec","acsc",	
                                "sqrt","rand","eye","zeros","ones","log10","log","exp"};


void trimend(string & x)
{
	for (int i = x.length() - 1; i<0; i--)
	{
		if (x[i] == ' ')
		{
			x.erase(i, 1);

			i--;
		}
		else
		{
			break;

		}
	}
}

void trimSpace(string & x)
{
	for (int i = 0; i<x.length(); i++)

	{
		if (x[i] == ' ')
		{
			x.erase(i, 1);

			i--;
		}

	}



}
void trimline(string & x)
{
	for (int i = 0; i<x.length(); i++)

	{
		if (x[i] == '\n' || (int)x[i]==13)
		{
			x.erase(i, 1);

			i--;
		}

	}



}


/*  get varible index from  the array   */
int index(string varName, const string * varibleNames, int NVar)
{

	for (int i = 0; i < NVar; i++)
	{


		if (varName == varibleNames[i])
		{



			return i;
		}

	}
	return (-1);

}

// Function to get weight of an operator. An operator with higher weight will have higher precedence.
int GetOperatorWeight(char op)
{
	int weight = -1;
	switch(op)
	{
	case '(':
        {weight = 0;
        break;}
    case '+':
	case '-':
		{weight = 1;

		break;}
	case '*':
	case '/':
    case '%':
		{weight = 2;
		break;}
	case '^':
	case '\'':
		{weight = 3;
		break;}
	}
	return weight;
}

// Function to perform an operation and return output.
int HasHigherPrecedence( string op1,  string op2)
{
	int op1Weight = GetOperatorWeight(op1[0]);
	int op2Weight = GetOperatorWeight(op2[0]);

	return op1Weight > op2Weight ?  1: 0;
}
int Isoperator(char a)
{
	if(a=='+'||a=='-'||a=='*'||a=='/'||a=='%'||a=='^'||a=='\'')return 1;
	else return 0;
}
CMatrix postfixCalculation(string  s)
{
		trimSpace(s);
    
	if(s[0]=='-'||s[0]=='+')
      {
    		s.insert(0,1,'0');

      }

    for(int i=1;i<s.length()-1;i++)
    {
    	if(s[i]=='\'')
    		{
    			 if(!(Isoperator(s[i+1])||s[i+1]==')'))
    			 	   s.insert(i+1,1,'*');

    		}
    	if((s[i]=='-'||s[i]=='+')&&s[i-1]=='(')
    		{
    			s.insert(i,1,'0');	
    		}
    }


    if(s.length()>1)
    {
    for(int i=1;i<(s.length()-2);i++)
       {
    	
    	if(s[i]==')'&&!(Isoperator(s[i+1])||s[i+1]=='\'')&&s[i+1]!=')'&&(s[i+1]!='.'&&Isoperator(s[i+2])))
    	{
    	  	s.insert(i+1,1,'*');
    	}

    	if(s[i]=='('&&(!Isoperator(s[i-1]))&&s[i-1]!='(')
    	{
    	  	s.insert(i,1,'*');
    	}
       }
    }

   


    

    int flag=1;
    int flagpower=0;
    stack <string> ope,post,ex;
    stack <CMatrix> result;
    int sum=0;
    char* infix = new char[s.length() + 1];
	strcpy(infix, s.c_str());
	char* spearators = "+/*-^% \'";
	char* token = strtok(infix, spearators);
	while (token)
	{ int counter=0;
	    sum+=strlen(token);
	    if(token[strlen(token)-1]=='.'){token[strlen(token)-1]=NULL;flagpower=1;}
	    while(token[0]=='('){ope.push("(");token++;}
	    while(token[strlen(token)-1]==')')
            {
                token[strlen(token)-1]=NULL;   counter++; 
            }

        post.push(token);
        while(counter)
           {
                while(ope.top()!="(")
                    {
                            post.push(ope.top());
                            ope.pop();
                    }
                    ope.pop();
                    counter--;
            }



         for(int i=0;i<1;i++)
            {
                if(ope.empty()) continue;
                 while((!HasHigherPrecedence(s.substr(sum,1),ope.top()))&&(!ope.empty()))
                        {
                            post.push(ope.top());
                            ope.pop();
                            if(ope.empty()) break;
                        }

            }

	    if(sum<s.length()) 
	    	{
	    		if((flagpower==1)&&(Isoperator(s[sum]))) ope.push(s.substr(sum,1)+".");
	    		else ope.push(s.substr(sum,1));
            }
        sum++;
		token = strtok(NULL, spearators);
		flag++;
		flagpower=0;
	}

	

while(!post.empty())
{
	ex.push(post.top());
	post.pop();
}


int L=0;
while(!ex.empty())
{
	string top=ex.top();
	int x=index(top,tempvaribleNames,tempNVar);
	int y=index(top,varibleNames,NVar);
	if(Isoperator(top[0]))
	{
		CMatrix A =result.top();
		result.pop();
        if(top[0]=='+')
        {
        	if(A.getn()==1) result.top()+=A(0,0);
        	else if(result.top().getn()==1) result.top()=A+result.top()(0,0);
        	else	result.top()+=A;
        }
        else if(top[0]=='-')
        {

        	if(A.getn()==1) result.top()-=A(0,0);
        	else if(result.top().getn()==1) result.top()=-A+result.top()(0,0);
        	else	result.top()-=A;
        }
        else if(top[0]=='*')
        {

        	if(A.getn()==1) result.top()=result.top()*A(0,0);
        	else if(result.top().getn()==1) result.top()=A*result.top()(0,0);
        	else	{
        		              if(top.length()==1)result.top()=result.top()*A;
        		              else
        		              {
        		              	if(A.getnR()!=result.top().getnR() ||A.getnC()!=result.top().getnC()) throw("error  Invalid Matrices Dim");
        		              	else
        		              	{
        		              		for(int r=0;r<A.getnR();r++) for(int c=0;c<A.getnC();c++)
        		              		 { result.top().setElement(r,c,result.top()(r,c)*A(r,c));}
        		              	}



        		              }
        		    }
        }
        else if(top[0]=='/')
        {
        	
        	if(A.getn()==1){ result.top()=result.top()*(1.0/A(0,0)); if(A(0,0)==0) L=1;} 
        	else if(result.top().getn()==1) {result.top()=A.elementDiv(result.top()(0,0)); if (result.top()(0,0)==0)L=1;}
        	else	{
        		              if(top.length()==1)result.top()=result.top()/A;
        		              else
        		              {
        		              	if(A.getnR()!=result.top().getnR() ||A.getnC()!=result.top().getnC()) throw("error  Invalid Matrices Dim");
        		              	else
        		              	{
        		              		for(int r=0;r<A.getnR();r++) for(int c=0;c<A.getnC();c++)
        		              		 result.top().setElement(r,c,result.top()(r,c)/A(r,c));
        		              	}



        		              }
        		    }

        }
        else if(top[0]=='%')
        {
        	int LS,RS;
        	LS=result.top()(0,0);
        	RS=A(0,0);
        	LS=LS%RS;
        	CMatrix H(LS);
        	result.top()=H;
        }

        else if(top[0]=='^')
        {    if (A.getn()==1)
        	{
        	if(top.length()>1) result.top()=result.top().powElement(A(0,0));
        	else if(result.top().getn()==1) result.top()=result.top().powElement(A(0,0));
        	else result.top()=result.top().pow_matrix(A(0,0));
       	    }
       	    else{
                    if(A.getnR()!=result.top().getnR() ||A.getnC()!=result.top().getnC()) throw("error:  Invalid Matrices Dim");
                    else
                    {

                    		for(int r=0;r<A.getnR();r++) for(int c=0;c<A.getnC();c++)
        		              		 result.top().setElement(r,c,pow(result.top()(r,c),A(r,c)));

                    }
       	    }

        }
        else if(top[0]=='\'')
        {
        	result.push(A.getTranspose());
        }
	}
    else if(y!=-1)
    {
    	result.push(Matrices[y]);
    }
    else if(x!=-1)
    {
        result.push(tempMatrices[x]);

    }

    else
    {
    	if(isdigit(top[0]) )
    		{
    			CMatrix m(atof(top.c_str()));
    	        result.push(m);
    	    }

    	else if(top.length()==1)
    	    {
                   throw("undefine varible");
    	    }
    	else
    	{
    		if(isdigit(top[1])&&top[0]=='.')
    		{
    			
    			CMatrix m(atof(top.c_str()));
    	        result.push(m);
    	        
    	    
    		}
    		else throw("undefine varible");
    	}
    }




    ex.pop();

}
if(L==1)cout<<"warning: you want to divide to zero "<<endl;

return result.top();



}






CMatrix functionExe(int i,CMatrix &m,int r=0,int c=0)
{
	switch(i)
	{   case 0:
	           {
	    	     return sinh(m);
	             break;

	    	   }
	         
	    case 1: 
	            {
	            	return cosh(m);
	                break;
	            }
	         
	     case 2:
	           {
	    	     return tanh(m);
	             break;

	    	   }
	    case 3:
	           {
	    	     return coth(m);
	             break;

	    	   }
	    case 4:
	           {
	    	     return sech(m);
	             break;

	    	   }
	    case 5:
	           {
	    	     return csch(m);
	             break;

	    	   }
	         
	    case 6: 
	            {
	            	return sin(m);
	                break;
	            }
	         
	    case 7:
	           {
	    	     return cos(m);
	             break;

	    	   }
	    case 8:
	           {
	    	     return tan(m);
	             break;

	    	   }
	    case 9:
	           {
	    	     return cot(m);
	             break;

	    	   }
	    case 10:
	           {
	    	     return sec(m);
	             break;

	    	   }
	         
	    case 11: 
	            {
	            	return csc(m);
	                break;
	            }
	         
	    case 12:
	            {
	    	     return asin(m);
	             break;

	    	    }
	    case 13:
	           {
	    	     return acos(m);
	             break;

	    	   }
	    case 14:
	           {
	    	     return atan(m);
	             break;

	    	   }
	    case 15:
	           {
	    	     return acot(m);
	             break;

	    	   }
	         
	    case 16: 
	            {
	            	return asec(m);
	                break;
	            }
	         
	     case 17:
	           {
	    	     return acsc(m);
	             break;

	    	   }
	    case 18:
	           {
	    	     return sqrt(m);
	             break;

	    	   }
	   case 19:
	           {
	           	
	    	     return CMatrix (r,c,3,0);
	             break;

	    	   }
            					                  
        case 20:
	           { 
	    	     return CMatrix (r,c,2,0);
	             break;

	    	   }
	         
	    case 21: 
	            {
	            	return CMatrix (r,c,0,0);
	                break;
	            }
	         
	     case 22:
	           {
	    	     return CMatrix (r,c,1,0);
	             break;

	    	   }
	    case 23:
	           {
	    	     return log(m);
	             break;

	    	   }
	    case 24:
	           {
	    	     return ln(m);
	             break;

	    	   }
	     case 25:
	           {
	    	     return exp(m);
	             break;

	    	   }
	         
	                					           
 
	         
	}  		
	         
}
CMatrix Execution(string s )
{

    trimSpace(s);
	trimline(s);
	for (int x=0;x<26;x++)
	{
    		while(s.find(functionsname[x])!=-1)
    			{	
    				    string func;
    				    int index=s.find(functionsname[x]);
    				    int counter=0;
    				    int flag=0;
        				for(int z=index;z<s.length();z++)
        					{	
            					if(s[z]=='(')
                					{	counter++;
                    					flag=1;
                    				}
            					if(s[z]==')')
                					{	counter--;
                    					flag=1;
                    				}
            					if(counter==0&&flag==1)
            						{
                						func=s.substr(index+functionsname[x].length()+1,z-index-functionsname[x].length()-1);
                						string func1;
                						CMatrix e(0);
                						if(x==19||x==20||x==21||x==22)
                						{
                							func1=func.substr(func.find(",")+1);
                							func.erase(func.find(","));
                							e=Execution(func1);
                						}
                							
                						CMatrix c=Execution(func);
                					
                						c=functionExe(x,c,c(0,0),e(0,0));
    									
                						tempMatrices[tempNVar]=c;
                						////put the name of temp matrix 
                						
                						
                						
                						
                						if(index==0)
                						{
                							s.erase(index,z-index+1);
                							s=tempvaribleNames[tempNVar]+s;
                						}
                						else if(z==(s.length()-1))
                						{
                							s.erase(index,z-index+1);
                							s=s+tempvaribleNames[tempNVar];
                						}
                						else 
                						{
                								s=s.substr(0,index)+tempvaribleNames[tempNVar]+s.substr(z+1);
                						}
                							tempNVar++;
                							break;
            						}


        					}	
    			}




	}

  

	return postfixCalculation(s);

}
void errordetection(string &s)
{

	int counter =0;
	if(((   Isoperator(s[0]) && (s[0]!='+'&&s[0]!='-') || (Isoperator(s[s.length()-1]) &&s[s.length()-1]!='\''))))  throw("error:invalid expression");
	for(int i=0;i<s.length()-1;i++)
	{

		if(Isoperator(s[i+1])&&Isoperator(s[i]))
		{
			if(s[i+1]=='*'&&s[i]=='*')
			{
				s[i]='.';s[i+1]='^';
			}
			else if(s[i]==s[i+1])
			{

				throw("error:invalid expression");
			}
			else if(s[i+1]=='+'&& s[i]!='\'')
			{
                       s.erase(i+1,1);
			}
			else if(s[i+1]=='-' && s[i]!='^')

			{
				    /*int num =0;
				    s.replace(i+1,1,"0-1*");
				    
                     int j=i+5;
                     int q=s.length();
				    for( ; j<q ;j++)
				    	{
				    		
				    		if(Isoperator(s[j]) && num==0)  { char w =s[j]; s.replace(j,1,")"+w); break;}
				    	    else if (s[j]=='(') num++;
				    	    else if (s[j]==')') num--;
				    	    if(j==(s.length()-1)) s=s+")"; 
				    	}	
                     */
                       throw("error: please use () around negative value ");
			}
			else{


					throw("error:invalid expression");
				}
		}


		if(s[i]=='(')
		{
			counter++;
			if(Isoperator(s[i+1]) && ((s[i+1]!='+')&&(s[i+1]!='-'))) throw("error:invalid expression");
		}
		else if(s[i]==')')
		{
			counter--;
			if(Isoperator(s[i-1])) throw("error:invalid expression");
		}

	}
	if(s[s.length()-1]==')') counter--;
	if(counter!=0) throw("error:invalid expression");
}

void  operationHandling (string s)
{

	trimSpace(s);
	trimline(s);
	if(s.length()==0) return;
	int printFlag=1;
	if(s.find(";")!=-1){printFlag=0; s.erase(s.find(";"));}
    
    errordetection(s);
	
	

	if(s.find("=")==-1)
	{
			if(printFlag==0)return;
			if(index(s,varibleNames,NVar)==-1)cout<<"ans ="<<endl;
			else cout<<s<<" ="<<endl;
			CMatrix result;
			result=Execution(s);
            cout<<result;
            cout<<endl;
            return;
	}
	else
	{
		string destinations,calculations;
		int u=s.rfind("=");
		destinations=s.substr(0,u+1);
		calculations=s.substr(u+1);
		CMatrix result;
		result=Execution(calculations);
        int destinationsNum=1;
        u=destinations.find("=");
        int flagz=0;
        while(u!=-1)
        {
        	int i=index(destinations.substr(0,u),varibleNames,NVar);
        	if(i!=-1) 
        		{
        			Matrices[i]=result;
        			if(printFlag==1 &&flagz==0) 
        				{
        					cout<<varibleNames[i]<<" ="<<endl;
        					flagz++;
        				}
        		}
        	else
        	{
        		Matrices[NVar]=result;
        		varibleNames[NVar]=destinations.substr(0,u);
        		if(printFlag==1 &&flagz==0) 
        				{
        					cout<<varibleNames[NVar]<<" ="<<endl;
        					flagz++;
        				}
        		NVar++;
        	}
        	destinations.erase(0,u+1);
        	if(destinations.length()<=0) u=-1;
        	else u=destinations.find("=");
        }
        if(printFlag==1)
        {
        	cout<<result<<endl;
        }
		tempNVar=0;
	}
}
//////////////////////////////////////////////////medhat amira aya/////////////////////////////////////////////////////////////////////////////////////


string modify(string& s);
void parse(string &s);
int findOPs(string& s, int pos);
int findSpaceBefore(string& s, int pos);
int findSpaceAfter(string& s, int pos);
string getElement(string s1);
void updateString(string& s1);
void trim(string& text);
int myfind(string& s, int pos);
void advancedTrim(string& s);
string advancedConcatination(string& s);
void trimbegin(string & x);






int main(int argc, char*argv[])
{

	
	 
	
	
				string op;

				if (argc == 1)
				{
					
					while (1)
					{
						try{
								cout<<">>";
								getline(cin, op);
								int m;
								m = op.find("[");
								if (m!=-1)
								{
									int h=0;
										for (int i=0;i<(op.length()) ;i++)
									    {
										   if(op[i]=='[') {h++;}
										   if(op[i]==']') {h--;}
									    }


									while (1)
									{
										
										if (h!=0)
										{
											string temp;
											getline(cin, temp);
											if(temp.length()!=0)
											{
												if(op[op.length()-1]=='\n') op = op + temp;
											    else op=op+"\r\n"+temp;
											
											    for (int i=0;i<(temp.length()) ;i++)
									                {
										                    if(temp[i]=='[') h++;
										                    if(temp[i]==']') h--;
									                }
									        }

										}
										else
										{
											break;
										}


									}

									parse(op);
										
								}
                                    
							    else  operationHandling(op);

								


							}
							catch( const char * error)
								{
									cout<< error<<endl;
								}

					}

				}
				else
				{
					ifstream infile(argv[1]);
					while (!infile.eof())
					{

						try 
						{
							getline(infile, op);
						    int m;
						    m = op.find("[");
						    if (m!=-1)
						       {
						       	    int h=0;
						       	    for (int i=0;i<(op.length()) ;i++)
									    {
										   if(op[i]=='[') h++;
										   if(op[i]==']')h--;
									    }

							        while (1)
							         {
							         	
								        
								        if (h!=0)
								          {
									         string temp;
									         getline(infile, temp);
									         if(temp.length()!=0)
											{
												if(op[op.length()-1]=='\n' ) op = op + temp;
											    else op=op+"\r\n"+temp;
											
											    for (int i=0;i<(temp.length()) ;i++)
									                {
										                    if(temp[i]=='[') h++;
										                    if(temp[i]==']') h--;
									                }
									        }
									         
								          }
								        else break;								
							         }
							         cout<<op<<endl;
                                     for(int i=0;i<op.length();i++) if((int)op[i]==13){op.erase(i,1);i--;}
							         parse(op); 
						        }

						    else operationHandling(op);


						}
						catch( const char * error)
						{
									cout<< error<<endl;
						}

						


					}
				}
				
			

			
		

	delete[] varibleNames;
	delete[] Matrices;
	delete[] tempMatrices;







}





//////////////////////////////////////////////////////////Matrix parsing/////////////////////////  
  
  void advancedTrim(string& s)
{
	for (int i = 0; i < s.length(); i++)
	{
		if (s[i] == ' ' && s[i + 1] == ' ')
		{
			s.erase(i, 1);
			i--;
		}
	}
}

void updateString(string& s1)                /*replace the part to send to alaa */
{
	string s;
	trim(s1);
	int pos = findOPs(s1, 0), pos1, pos2;
	while (pos != -1)
	{
		pos1 = findSpaceBefore(s1, pos);
		pos2 = findSpaceAfter(s1, pos);
		s = s1.substr(pos1 + 1, pos2 - pos1 - 1);
		CMatrix g=Execution(s);
		if(g.getn()==1) {s=g.sendString(); s.erase(s.find("["),1);s.erase(s.find("]"),1);}
		else s=s=g.sendString();
		s1 = s1.replace(pos1 + 1, pos2 - pos1 - 1, s); 
		pos = findOPs(s1, pos+1);
	}
	
}

void trim(string& text)      //////////trimming spaces around operations//////////
{
	int n, startpos = 0;
	int pos;
	pos = myfind(text, startpos);
	do
	{
		n = 0;
		if (pos != -1)
		{
			if (text[pos + 1] == ' ')
			{
				text.erase(pos + 1, 1);
				startpos = pos + 1;
				n++;
			}
			if (text[pos - 1] == ' ')
			{
				text.erase(pos - 1, 1);
				startpos = pos;
				n++;
			}
			if (n == 0)
				startpos = pos + 1;
			pos = myfind(text, startpos);
		}
	} while (pos != -1);
}

int myfind(string& s, int pos) ////////// find the operations for trim //////////////////////
{
	for (int i = pos; i < s.length(); i++)
	{

		if (s[i] == '+' || s[i] == '-' || s[i] == '*' || s[i] == '/' || s[i] == '^' || s[i] == ';' || s[i] == ','
			|| (s[i] == '.' && s[i + 1] == '+') || (s[i] == '.' && s[i + 1] == '-') || (s[i] == '.' && s[i + 1] == '*')
			|| (s[i] == '.' && s[i + 1] == '/') || (s[i] == '.' && s[i + 1] == '^'))
			return i;
	}

	return -1;
}
int findOPs(string& s, int pos)         // find operations to send to be calculated 
{
	for (int i = pos; i < s.length(); i++)
	if (s[i] == '+' || s[i] == '-' || s[i] == '*' || s[i] == '/' || s[i] == '^' || (s[i] == 's' && s[i + 1] == 'i' && s[i + 2] == 'n')
		|| (s[i] == 'c' && s[i + 1] == 'o' && s[i + 2] == 's') || (s[i] == 't' && s[i + 1] == 'a' && s[i + 2] == 'n')
		|| (s[i] == 'l' && s[i + 1] == 'o' && s[i + 2] == 'g') || (s[i] == 'l' && s[i + 1] == 'n')
		|| (s[i] == 's' && s[i + 1] == 'q' && s[i + 2] == 'r' && s[i + 3] == 't'))
		return i;
	return -1;
}
int findSpaceBefore(string& s, int pos)         // find begining of operation
{
	for (int i = pos; i >= 0; i--)
	if (s[i] == ' ' || s[i] == '[' || s[i] == ']' || s[i] == ';')
		return i;
	return -1;
}
int findSpaceAfter(string& s, int pos)         // find the end of operation
{
	for (int i = pos; i < s.length(); i++)
	if (s[i] == ' ' || s[i] == '[' || s[i] == ']' || s[i] == ';')
		return i;
	return -1;
}
void parse(string &operation)             /////take the string and return it in the form of phase 1
{
	string varName, s2, s;
	int pos, flag = 0;
	CMatrix x;
	if (operation.find("=") != -1)
	{
		varName = operation.substr(0, operation.find("="));
		trimSpace(varName);
		s2 = operation.substr(operation.find("=") + 1);
		if (s2.find("\r\n") != -1) s2 = s2.replace(s2.find("\r\n"), 2, ";");

		if (s2[s2.length() - 1] == ';')   // to know whether to print or not
		{
			flag = 1;
			s2.erase(s2.length() - 1, 1);
		}

		trimend(s2);					//removing extra spaces and update string
		trimbegin(s2);
		advancedTrim(s2);
		updateString(s2);


		for (int i = 0; i < NVar; i++)   //search if there's an existant matrix variable name (A or B etc.) in the string 
		{
			pos = s2.find(varibleNames[i]);
			while (pos != -1)
			{
				s2.replace(pos, varName.length(), Matrices[i].sendString());
				pos = s2.find(varibleNames[i]);
			}
		}
        
		s2 = advancedConcatination(s2);   // the string is in the standard format of matrix
		
		x.copy(s2);
         
		if(index(varName,varibleNames,NVar)!=-1)
		{
			Matrices[index(varName,varibleNames,NVar)]=x;
		}
		else 
		{
             Matrices[NVar]=x;
             varibleNames[NVar]=varName;
             NVar++;
		}
		if (flag == 0){ cout << varName<<" ="<<endl<<x<<endl; } /*print the matrix when there is no semicolon in the end of the operation*/
	}

	//if just the Matrix name is written to be printed
	else
	{
		varName = operation;
		trimSpace(varName);
		trimend(varName);
		for (int i = 0; i <= NVar; i++)
		{
			if (varibleNames[i] == varName)
				cout << Matrices[i].getString() << endl;
			else throw("Error:you try to print undefined Matrix");
		}
	}

	
}

string modify(string& s)
{
	CMatrix A;
	CMatrix B;
	string s1, s2, s3;
	int i;
	int j;
	int k;

	while (s.find("] [") != -1 || s.find("],[") != -1)
	{
		i = s.find("] [");
		k = s.find("],[");
		if (i != -1)
		{
			s1 = s.substr(i + 2);
		}
		else if (k != -1)
		{
			s1 = s.substr(k + 2);
		}
		j = s1.find("]");
		s1 = s1.substr(0, j + 1);
		A.copy(s1);

		if (i != -1)
		{
			s2 = s.substr(0, i + 1);
		}
		else if (k != -1)
		{
			s2 = s.substr(0, k + 1);
     }
		j = s2.rfind("[");
		s2 = s2.substr(j);
		B.copy(s2);
		B.concatinate(A);
		s = s.replace(j, s1.length() + s2.length() + 1, B.sendString());
	}
	while (s.find("[") != -1)
	{
		s.erase(s.find("["), 1);
	}
	while (s.find("]") != -1)
	{
		s.erase(s.find("]"), 1);
	}
	s = s + "]";
	s = '[' + s;
	if (s[1] == ' ') s.erase(1, 1);
	return s;

}

string advancedConcatination(string &s)  // s= [[1.2 2.3; 1 2.3; [1.3 2.4;4.6 1.3]],[3.2 1;-7.8 2;-3.2 3; 1.2 4]]
{
	int i, j, k, l;
	string s1, s2;
	CMatrix A;
	i = s.find("] [");
	if (i != -1) s = s.replace(i + 1, 1, ",");
	i = s.find("],[");
	if (i < 0)
	{
		A.copy(s); return A.sendString();
	}
	while (i != -1)
	{
		s1 = s.substr(0, i);						// s1= [[1.2 2.3; 1 2.3; [1.3 2.4;4.6 1.3]
		j = s1.rfind('[');
		//if (j < 0) throw ("Error: invalid Matrix");
		k = s1.find(']', j);
		while (j != -1 && k != -1)
		{
			s1.erase(k, 1);
			s1.erase(j, 1);							//s1= [[1.2 2.3; 1 2.3; 1.3 2.4;4.6 1.3
			j = s1.rfind('[');
			k = s1.find(']', j);
		}
		s.replace(0, i, s1);						// s= [[1.2 2.3; 1 2.3; 1.3 2.4;4.6 1.3],[3.2 1;-7.8 2;-3.2 3; 1.2 4]]
		i = s.find("],[", i + 2);
	}
	i = s.rfind("],[");
	s1 = s.substr(i + 3);
	l = s1.length();
	j = s1.rfind('[');
	k = s1.find(']', j);
	while (j != -1 && k != -1)
	{
		s1.erase(k, 1);
		s1.erase(j, 1);							//s1= [[1.2 2.3; 1 2.3; 1.3 2.4;4.6 1.3
		j = s1.rfind('[');
		k = s1.find(']', j);
	}
	s.replace(i + 3, l+3, s1);

	s2 = modify(s);
	A.copy(s2);
	return A.sendString();
}


  void trimbegin(string & x)
{
	for (int i = 0; i<x.length() - 1; i++)
	{
		if (x[i] == ' ')
		{
			x.erase(i, 1);
			i--;
		}
		else
			break;
	}
}
  
  
  /////////////////////////////////end of Matrix Parsing ////////////////////////////////////////////
