#include <iostream>
#include <cmath>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <sstream>
using namespace std;

struct polarNum;
struct currentPath;
class connection;
class node;

polarNum polarParse(float,float);

struct polarNum 
{
	float l;
	float a;
	float real()
	{
		return abs(l*cos(a))<abs(l*sin(a))*1E-5?0:l*cos(a);
	}
	float imagenary()
	{
		return abs(l*sin(a))<abs(l*cos(a))*1E-5?0:l*sin(a);
	}
	polarNum operator+(int b)
	{
		polarNum temp;
		temp.l = this->l+b;
		temp.a = this->a;
		return temp;
	}
	polarNum operator+(polarNum b)
	{
		return polarParse(this->l*cos(this->a)+b.l*cos(b.a),this->l*sin(this->a)+b.l*sin(b.a));
	}
	polarNum operator-(polarNum b)
	{
		return polarParse(this->l*cos(this->a)-b.l*cos(b.a),this->l*sin(this->a)-b.l*sin(b.a));
	}
	polarNum operator/(polarNum b)
	{
		polarNum temp;
		temp.l = this->l/b.l;
		temp.a = this->a-b.a;
		return temp;
	}
	polarNum operator*(polarNum b)
	{
		polarNum temp;
		temp.l = this->l*b.l;
		temp.a = this->a+b.a;
		return temp;
	}
	polarNum operator/(float b)
	{
		polarNum temp;
		temp.l = this->l/b;
		temp.a = this->a;
		return temp;
	}
	polarNum operator*(float b)
	{
		polarNum temp;
		if(b<0)
		{
			temp.l = this->l*b*(-1);
			temp.a = this->a + M_PI;
		}
		else
		{
			temp.l = this->l*b;
			temp.a = this->a;
		}
		return temp;
	}
	polarNum operator=(double b)
	{
		this->l = b;
		this->a = 0;
		return *this;
	}
	bool operator==(polarNum b)
	{
		if(this->l==b.l && this->a == b.a)
			return true;
		return false;
	}
	bool operator!=(polarNum b)
	{
		if(this->l==b.l && this->a == b.a)
			return false;
		return true;
	}
	bool operator==(double b)
	{
		if(this->l==b)
			return true;
		return false;
	}
	bool operator!=(double b)
	{
		if(this->l==b)
			return false;
		return true;
	}
};

class node
{
private:
	int connectionCount;
	int connectionCounter;
	int selfNum;
	connection **connectionArray;

public:
	static node* nodes;
	static int count;
	node();
	~node();
	void setConnectionCount(int,int);
	void addConnectionDetails(connection*);
	connection* getNotTraversedConnection();
	connection* getConnectionTo(node*,connection* const);
	connection* getNextConnection();
	int getSelfNum();
	void resetConnectionCounter();
};

node* node::nodes = NULL;
int node::count = 0;

class connection
{
private:
	node *fromNode;
	node *toNode;
	bool traversed;
	int *variable;
	polarNum impedence;
	polarNum voltage;
	polarNum current;
public:
	static connection* connections;
	static int count;
	connection();
	~connection();
	void setDetails(node*,node* ,polarNum,polarNum,polarNum);
	void printData();
	void setTraversed(node*,int);
	bool isTraversed();
	void setPotentialEquation(polarNum[],int);
	void setCurrentEquation(polarNum[]);
	node* getOtherNode(node*);
	int getFromNodeNumber()
	{
		return fromNode->getSelfNum();
	}
	int getToNodeNumber()
	{
		return toNode->getSelfNum();
	}
	bool isCurrentSourcePresent()
	{
		if(current==0) return false;
		return true;
	}
	int getCurrentDirectionFor(int i)
	{
		return variable[i];
	}
};

connection* connection::connections = NULL;
int connection::count = 0;

struct currentPath
{
	node *through;
	connection *from;
	connection *to;
	currentPath *next;
	currentPath *supposeNext;
	int selfNum;
	bool traversed;
	static int pathCountActual;
	static int maxLoops;
	static int pathCountSuppose;
	static int currentCounter;
	static currentPath **allPaths;
	currentPath(node*,currentPath*);
	static currentPath* make(node*,currentPath*);
	static currentPath* start;
	static currentPath* startNewTraversePath(node*);
	bool getNextNode();
	void solidifyPath();
	void convertSupposeToActual();
	static void startGettingNewEquation(polarNum[]);
	bool traverseToGetEquation(polarNum[],int);
	static bool onGettingCurrentSource(polarNum[],currentPath*,int);
	bool hasOccured(node*);
	void setTraversedFalse()
	{
		traversed = false;
	}
	bool isTraversed()
	{
		return traversed;
	}
	void setTraversedTrue()
	{
		traversed = true;
	}
	~currentPath();
	static void printActualCurrentPath();
	static void printSupposeCurrentPath();
};

class matrixSolver
{
private:
	polarNum **matA;
	int variableCount;
	int matrixCounter;
public:
	matrixSolver(int nVar)
	{
		variableCount = nVar;
		matrixCounter = 0;
		matA = new polarNum*[nVar];
		for (int i = 0; i < nVar; ++i)
		{
			matA[i] = new polarNum[nVar+1];
		}
	}
	void printMatrix()
	{
		cout<<"\nprinting matrix\n";
		for (int i = 0; i < variableCount; ++i)
		{
			for (int j = 0; j < variableCount+1; ++j)
			{
				cout<<matA[i][j].l<<" ";
			}
			cout<<endl;
		}
	}
	void addEquations(double *row)
	{
		for (int i = 0; i < variableCount+1; ++i)
		{
			matA[matrixCounter][i] = row[i];
		}
		matrixCounter++;
	}
	void addEquations(polarNum *row)
	{
		for (int i = 0; i < variableCount+1; ++i)
		{
			matA[matrixCounter][i] = row[i];
		}
		matrixCounter++;
	}
	void solve()
	{
		for (int i = 0; i < variableCount; ++i)
		{
			if(matA[i][i] == 0)
			{
				int j;
				for (j = i; j < variableCount; ++j)
				{
					if(matA[j][i]!=0)
						break;
				}
				if(j==variableCount && matA[i][variableCount]==0)
				{
					cout<<"Infinite Solutions are possible\n\n";
					break;
				}
				if(j==variableCount && matA[i][variableCount]!=0)
				{
					cout<<"No solution is possible\n\n";
					break;
				}
				polarNum *temp = matA[j];
				matA[j] = matA[i];
				matA[i] = temp;
			}
			for (int k = variableCount; k >= i; --k)
			{
				matA[i][k] = matA[i][k] / matA[i][i];
			}
			for (int j = 0; j < variableCount; ++j)
			{
				if(j == i)
						continue;
				for (int l = variableCount; l >=i; --l)
				{
					matA[j][l] = matA[j][l] - matA[i][l]*matA[j][i];
				}
			}
		}
	}
	polarNum getSolutionFor(int _var)
	{
		return matA[_var][variableCount];
	}
};

node::node()
{
	connectionCounter = 0;
}

void node::setConnectionCount(int n,int a)
{
	this->connectionArray = new connection*[n];
	this->connectionCount = n;
	this->selfNum = a;
}

void node::addConnectionDetails(connection *_new)
{
	connectionArray[connectionCounter] = _new;
	connectionCounter++;
	if(connectionCounter == connectionCount)
		connectionCounter = 0;
}

void node::resetConnectionCounter()
{
	connectionCounter = 0;
}

connection* node::getNotTraversedConnection()
{
	for (int i = 0; i < connectionCount; ++i)
	{
		if (!connectionArray[i]->isTraversed())
		{
			return connectionArray[i];
		}
	}
	return NULL;
}

connection* node::getConnectionTo(node* _node,connection* const notThis = NULL)
{
	for (int i = 0; i < connectionCount; ++i)
	{
		if(connectionArray[i] != notThis && connectionArray[i]->getOtherNode(this) == _node)
			return connectionArray[i];
	}
	return NULL;
}

connection* node::getNextConnection()
{
	if(connectionCounter==connectionCount)
		return NULL;
	return connectionArray[connectionCounter++];
}

int node::getSelfNum()
{
	return this->selfNum;
}

connection::connection()
{
	traversed = false;
	variable = NULL;
}

void connection::setDetails(node *prev, node *next,polarNum imp,polarNum volt,polarNum cur)
{
	this->fromNode = prev;
	this->toNode = next;
	this->impedence = imp;
	this->voltage = volt;
	this->current = cur;
}

void connection::printData()
{
	cout<<impedence.l<<"  "<<impedence.a/M_PI*180<<" : "<<voltage.l<<"  "<<voltage.a/M_PI*180<<" : "<<current.l<<"  "<<current.a/M_PI*180;
}

void connection::setTraversed(node *_from,int currentNumber)
{
	this->traversed=true;
	if(variable==NULL){
		variable = new int[currentPath::maxLoops];
		for (int i = 0; i < currentPath::maxLoops; ++i)	variable[i]=0;
	}
	variable[currentNumber] = _from==fromNode?1:-1;
}

void connection::setPotentialEquation(polarNum equationArray[],int currentNumber)
{
	int direct;
	if(currentNumber >= 0)
		if(variable[currentNumber] == -1)	direct = -1;
		else	direct = 1;
	else{
		if(variable[currentNumber * -1] == -1)	direct = 1;
		else	direct = -1;	
	}

	for (int i = 0; i < currentPath::maxLoops; ++i)
		equationArray[i]=equationArray[i]+impedence*variable[i]*float(direct);
	equationArray[currentPath::maxLoops] = equationArray[currentPath::maxLoops]+voltage*float(direct)*-1;
}

void connection::setCurrentEquation(polarNum equationArray[])
{
	for (int i = 0; i < currentPath::maxLoops; ++i)
		equationArray[i]=equationArray[i]+variable[i];
	equationArray[currentPath::maxLoops] = current;
}

bool connection::isTraversed()
{
	return this->traversed;
}

node* connection::getOtherNode(node* _node)
{
	return _node==toNode?fromNode:toNode;
}

currentPath* currentPath::start = NULL;
int currentPath::pathCountActual = 0;
int currentPath::pathCountSuppose = 0;
int currentPath::currentCounter = 0;
int currentPath::maxLoops = 0;
currentPath **currentPath::allPaths = NULL;

currentPath::~currentPath()
{
	if(this->next!=NULL)
		delete this->next;
}

currentPath::currentPath(node* _through,currentPath* _prev)
{
	this->through = _through;
	this->supposeNext = NULL;
	this->next = NULL;
	selfNum = currentCounter;
	traversed = false;
	if(_prev == NULL)
	{
		start = this;
		pathCountActual = 0;
		pathCountSuppose = 1;
		this->from = NULL;
	}
	else
	{
		_prev->supposeNext = this;
		this->from = _prev->to;
	}
}

currentPath* currentPath::make(node* _through,currentPath* _prev = NULL)
{
	return new currentPath(_through,_prev);
}

bool currentPath::hasOccured(node* _node)
{
	currentPath *temp = start;
	while(temp!=this){
		if(temp->through == _node)
			return true;
		temp = temp->supposeNext;
	}
	return false;
}

void currentPath::printActualCurrentPath()
{
	currentPath* temp = start;
	cout<<"current Path : \n";
	while(temp!=NULL)
	{
		cout<<temp->through->getSelfNum()<<" > ";
		temp = temp->next;
	}
	cout<<endl<<endl;
}

void currentPath::printSupposeCurrentPath()
{
	currentPath* temp = start;
	while(temp!=NULL)
	{
		temp = temp->supposeNext;
	}
	cout<<endl<<endl;
}

void currentPath::convertSupposeToActual()
{
	currentPath* temp=start;
	while(temp->supposeNext == temp->next){
		temp = temp->next;
	}
	if(temp->next!=NULL)
		delete temp->next;
	temp->next = temp->supposeNext;
	while(temp->supposeNext!=NULL){
		temp->next = temp->supposeNext;
		temp = temp->next;
	}
	pathCountActual = pathCountSuppose;
}

void currentPath::solidifyPath()
{
	currentPath* temp = start->next;
	start->to = start->through->getNotTraversedConnection();
	start->to->setTraversed(start->through,selfNum);
	while(temp->next!=NULL)
	{
		temp->to = temp->through->getConnectionTo(temp->next->through);
		temp->to->setTraversed(temp->through,selfNum);
		temp = temp->next;
	}
	temp->to = temp->through->getConnectionTo(start->through,start->to);
	temp->to->setTraversed(temp->through,selfNum);
	temp->next = start;
	start->from = temp->to;
}

bool currentPath::getNextNode()
{
	connection *temp;
	if( pathCountSuppose > pathCountActual && pathCountActual != 0)
	{
		pathCountSuppose--;
		return false;
	}
	if(this==start)
	{
		pathCountSuppose++;
		if(currentPath::make(this->through->getNotTraversedConnection()->getOtherNode(this->through),this)->getNextNode())
			return true;
		else
			return false;
	}
	if((temp=this->through->getConnectionTo(start->through,start->through->getNotTraversedConnection()))!=NULL)
	{
		convertSupposeToActual();
		pathCountSuppose--;
		return true;
	}
	this->through->resetConnectionCounter();
	while((temp=this->through->getNextConnection())!=NULL)
	{
		if(!hasOccured(temp->getOtherNode(this->through)))
		{
			pathCountSuppose++;
			if(!currentPath::make(temp->getOtherNode(this->through),this)->getNextNode())
			{
				delete this->supposeNext;
				this->supposeNext = NULL;
			}
		}
	}
	pathCountSuppose--;
	return true;
}


currentPath* currentPath::startNewTraversePath(node* _start)
{
	if(_start->getNotTraversedConnection()==NULL)
		return NULL;
	if(currentPath::make(_start)->getNextNode())
	{
		start->solidifyPath();
		return start;
	}
	else
	{
		delete start;
		return NULL;
	}
}

bool currentPath::onGettingCurrentSource(polarNum arr[],currentPath *path,int direction)
{
	int j,flag;
	for (j = 0; j < currentPath::maxLoops; ++j)
	{
		currentPath * temp = currentPath::allPaths[j];
		do
		{
			flag = 0;
			if(temp->to == path->to && temp != path)
			{
				temp->setTraversedTrue();
				if(temp->through == path->through)	temp->next->traverseToGetEquation(arr,direction*-1);
				else	temp->next->traverseToGetEquation(arr,direction);
				flag = 1;
				break;
			}
			temp = temp->next;
		}while(temp!=currentPath::allPaths[j]);
		if(flag == 1)	break;
	}
	if(flag==0)
	{
		currentPath *temp = path;
		do{
			temp->setTraversedTrue();
			temp = temp->next;
		}while(temp!=path);
		for(int i = 0 ;i<currentPath::maxLoops;i++)	arr[i].l = arr[i].a = 0;
		return false;
	}
	return true;
}

bool currentPath::traverseToGetEquation(polarNum arr[],int direction = 1)
{
	currentPath *temp = this;
	do
	{
		if(!temp->to->isCurrentSourcePresent())
		{
			if(direction == 1)	temp->to->setPotentialEquation(arr,this->selfNum);
			else	temp->to->setPotentialEquation(arr,this->selfNum*-1);
		}
		else
			if(!temp->isTraversed())
				if(!currentPath::onGettingCurrentSource(arr,temp,direction))
					return false;
		temp->setTraversedTrue();
		temp = temp->next;
	}while(temp!=this);
	return true;
}

void currentPath::startGettingNewEquation(polarNum arr[])
{
	int i;
	static int num = 0;
	for (i = 0; i < currentPath::maxLoops; ++i)
		if(!currentPath::allPaths[i]->isTraversed())
			break;
	if(i < currentPath::maxLoops)
	{
		if(!currentPath::allPaths[i]->traverseToGetEquation(arr))
		{
			cout<<"single current source found\n";
			currentPath::startGettingNewEquation(arr);
		}
	}
	else
	{
		cout<<"current eqution\n\n";
		for(i = num;i<connection::count;i++)
			if(connection::connections[i].isCurrentSourcePresent()){
				connection::connections[i].setCurrentEquation(arr);
				break;
			}
		num = i;
	}
}

polarNum polarParse(float a, float b)
{
	polarNum temp;
	temp.l = sqrt(a*a+b*b);
	if(a==0&&b!=0)
		temp.a = (b>0?1:-1)*M_PI/2;
	else if(a==0&&b==0)
		temp.a = 0;
	else
	{
		temp.a = atan(b/a);
		if(a<0)
		{
			temp.a+=M_PI;
		}
	}
	return temp;
}

vector<string> split(string str,char delim = ' ')
{
	vector<string> arr;
	vector<char> temp;
	//arr.push_back()
	for (int i = 0; i < str.length(); ++i)
	{
		if(str[i] != delim)
			temp.push_back(str[i]);
		else
		{
			if(temp.size()!=0)
			arr.push_back(string(temp.begin(),temp.end()));
			temp.clear();
		}
	}
	arr.push_back(string(temp.begin(),temp.end()));
	//cout<<arr.size()<<"  "<<arr.front()<<endl;
	return arr;
}

polarNum parseAwesomeData(string data)
{
	polarNum num;
	if(data.find("ang") != string::npos)
	{
		num.l = atoi(data.substr(0,data.find("ang")).c_str());
		num.a = atoi(data.substr(data.find("ang")+4,data.length()-1).c_str())*M_PI/180.0;
	}
	else
	{
		float re=0,im=0;
		string normalData;
		normalData = data;
		istringstream(normalData)>>re>>im;
		num = polarParse(re,im);
	}
	return num;
}

struct nodeData{int selfNum,connectionsCount;};
struct connectionData{int nodeFrom,nodeTo;polarNum volt,resis,curr;};

void setAwesomeData(vector<nodeData> nodes,vector<connectionData> connections,vector<int> nodeNumber)
{
	node::count = nodes.size();
	node::nodes = new node[node::count];
	connection::count=connections.size();
	connection::connections = new connection[connection::count];
	for(int i=0;i<node::count;i++)
		node::nodes[i].setConnectionCount(nodes[i].connectionsCount,nodes[i].selfNum);
	for (int i = 0; i < connection::count; ++i)
	{
		int nodeFrom = find(nodeNumber.begin(),nodeNumber.end(),connections[i].nodeFrom) - nodeNumber.begin();
		int nodeTo = find(nodeNumber.begin(),nodeNumber.end(),connections[i].nodeTo) - nodeNumber.begin();
		connection::connections[i].setDetails(&node::nodes[nodeFrom],&node::nodes[nodeTo],connections[i].resis,connections[i].volt,connections[i].curr);
		node::nodes[nodeFrom].addConnectionDetails(&connection::connections[i]);
		node::nodes[nodeTo].addConnectionDetails(&connection::connections[i]);
	}
	currentPath::maxLoops = connection::count - node::count + 1;
	currentPath::allPaths = new currentPath*[currentPath::maxLoops];
}

void getDataAwesome()
{
	vector<string> data,tempVec;
	vector<nodeData> nodes;
	vector<connectionData> connections;
	vector<int> nodeNumber;
	nodeData nodeTemp; connectionData connectionTemp;
	string temp;
	//cin.ignore();	//Because cin leaves delimiter character in the buffer and hence any next getline or anything takes in the delimiter character and fills the first space
	cout<<"Enter the data as node1 node2 and different params (avoid spaces)\n\
	Eg. 1 2 10+2jV 10ang(10)R 3+2jA \n\
	Type 'solve' to compute\n\n";
	do{
		cout<<"connection #"<<data.size()+1<<" : ";
		getline(cin,temp);
		if(temp.find("solve")==string::npos)
		{
			if(temp.find_first_not_of(' ') != string::npos) //checking for empty string
				data.push_back(temp);
		}
	}while(temp.find("solve")==string::npos);
	connection::count = data.size();
	for (vector<string>::iterator i = data.begin(); i != data.end(); ++i)
	{
		//cout<<*i<<endl;
		connectionTemp.volt = connectionTemp.curr = connectionTemp.resis = 0;
		tempVec = split(*i);
		vector<int>::iterator iteratorTemp;
		if((iteratorTemp = find(nodeNumber.begin(),nodeNumber.end(),atoi(tempVec[0].c_str())))==nodeNumber.end()){
			nodeNumber.push_back(atoi(tempVec[0].c_str()));
			nodeTemp.selfNum = atoi(tempVec[0].c_str());
			nodeTemp.connectionsCount = 1;
			nodes.push_back(nodeTemp);
		} //istringstream(tempVec[0])>>temp;
		else
			nodes[iteratorTemp - nodeNumber.begin()].connectionsCount++;
		if((iteratorTemp = find(nodeNumber.begin(),nodeNumber.end(),atoi(tempVec[1].c_str())))==nodeNumber.end()){
			nodeNumber.push_back(atoi(tempVec[1].c_str()));
			nodeTemp.selfNum = atoi(tempVec[1].c_str());
			nodeTemp.connectionsCount = 1;
			nodes.push_back(nodeTemp);
		}
		else
			nodes[iteratorTemp - nodeNumber.begin()].connectionsCount++;
		for(int j = 2;j<tempVec.size();j++)
		{
			switch(tempVec[j][tempVec[j].length()-1])
			{
				case 'v':
				case 'V':
						//cout<<parseAwesomeData(tempVec[j].substr(0,tempVec[j].length()-1)).real()<<"+"<<parseAwesomeData(tempVec[j].substr(0,tempVec[j].length()-1)).imagenary()<<endl;
						connectionTemp.volt = connectionTemp.volt + parseAwesomeData(tempVec[j].substr(0,tempVec[j].length()-1));
					break;
				case 'a':
				case 'A':
						connectionTemp.curr = connectionTemp.curr + parseAwesomeData(tempVec[j].substr(0,tempVec[j].length()-1));
					break;
				case 'r':
				case 'R':
						connectionTemp.resis = connectionTemp.resis + parseAwesomeData(tempVec[j].substr(0,tempVec[j].length()-1));
					break;
				default: cout<<"error yo "<<tempVec[j][tempVec[j].length()-1]<<endl;
			}	
		}
		connectionTemp.nodeFrom = atoi(tempVec[0].c_str());
		connectionTemp.nodeTo = atoi(tempVec[1].c_str());
		connections.push_back(connectionTemp);
		tempVec.clear();
	}
	cout<<endl;
	setAwesomeData(nodes,connections,nodeNumber);
	for(int i = 0; i < nodes.size(); ++i)
		cout<<"Node #"<<i+1<<" : "<<nodes[i].selfNum<<" "<<nodes[i].connectionsCount<<endl;
	cout<<endl;
	for (int i = 0; i < connections.size(); ++i)
		cout<<"Connection #"<<i+1<<" : "<<connections[i].nodeFrom<<" "<<connections[i].nodeTo<<" "<<connections[i].resis.real()<<"+"<<connections[i].resis.imagenary()<<"j "<<connections[i].volt.real()<<"+"<<connections[i].volt.imagenary()<<"j "<<connections[i].curr.real()<<"+"<<connections[i].curr.imagenary()<<"j "<<endl;
	cout<<endl;
}
/*
void getData()
{
	int temp;
	
	cout<<"Enter the number of nodes : ";
	cin>>node::count;
	node::nodes = new node[node::count];
	
	cout<<"Enter the number of connections per node\n";
	for(int i=0;i<node::count;i++)	//setting node deatils
	{
		cout<<"Node #"<<i+1<<" : ";
		cin>>temp;
		node::nodes[i].setConnectionCount(temp,i);
		connection::count+=temp;
	}
	connection::count/=2;
	connection::connections = new connection[connection::count];
	cout<<"Enter the connection details as 'nodeA nodeB Resistance Inductance Capacitance VoltageSourcePeak VoltageSourcePhase CurrentSourcePeak CurrentSourcePhase'\n";
	for(int i=0;i<connection::count;i++)	//setting connection::connections details
	{
		int nodeA,nodeB;
		float Resistance,Inductance,Capacitance,VoltageSourcePeak,VoltageSourcePhase,CurrentSourcePeak,CurrentSourcePhase;
		cout<<"Connection #"<<i+1<<" : ";
		cin>>nodeA>>nodeB>>Resistance>>Inductance>>Capacitance>>VoltageSourcePeak>>VoltageSourcePhase>>CurrentSourcePeak>>CurrentSourcePhase;
		connection::connections[i].setDetails(&node::nodes[nodeA-1],&node::nodes[nodeB-1],polarParse(Resistance,Inductance-Capacitance),polarNum({VoltageSourcePeak,VoltageSourcePhase/180*M_PI}),polarNum({CurrentSourcePeak,CurrentSourcePhase/180*M_PI}));
		node::nodes[nodeA-1].addConnectionDetails(&connection::connections[i]);
		node::nodes[nodeB-1].addConnectionDetails(&connection::connections[i]);
	}
	currentPath::maxLoops = connection::count - node::count + 1;
	currentPath::allPaths = new currentPath*[currentPath::maxLoops];
	
}*/

matrixSolver* evaluateData()
{
	currentPath *one,*two;
	for (int i = 0; i < node::count; ++i)
	{
		while((one = currentPath::startNewTraversePath(&node::nodes[i]))!=NULL)
		{
			currentPath::allPaths[currentPath::currentCounter++] = one;
			cout<<"Path #"<<currentPath::currentCounter<<" : ";
			two = one;
			do
			{
				cout<<two->through->getSelfNum()<<" > ";
				two = two->next;
				two->setTraversedFalse();
			} while (two!=one);
			cout<<two->through->getSelfNum()<<endl<<endl;
		}
	}
	polarNum arr[currentPath::maxLoops+1];
	matrixSolver *mat = new matrixSolver(currentPath::maxLoops);
	for (int j = 0; j < currentPath::maxLoops; ++j)
	{
		for (int i = 0; i < currentPath::maxLoops+1; ++i)
		{
			arr[i].l=0;
			arr[i].a=0;
		}
		currentPath::startGettingNewEquation(arr);
		mat->addEquations(arr);
		/*for (int k = 0; k <= currentPath::maxLoops; ++k)
		{
			cout<<arr[k].real()<<" + ";
		}
		cout<<endl;*/ 
	}
	mat->solve();
	return mat;
}

void displayResult(matrixSolver *mat)
{
	for (int i = 0; i < currentPath::maxLoops; ++i)
	{
		cout<<"Current in loop #"<<i+1<<" = "<<(mat->getSolutionFor(i)).real()<<" + "<<(mat->getSolutionFor(i)).imagenary()<<"j"<<endl;
	}
	cout<<"\n\n\n\nCurrent from one node to other are given below :-\n\n";
	for(int i = 0;i<connection::count;i++)
	{
		polarNum sum = {0,0};
		for(int j = 0;j<currentPath::maxLoops;j++)
		{
			sum = sum + mat->getSolutionFor(j) * connection::connections[i].getCurrentDirectionFor(j);
		}
		cout<<"Node #"<<connection::connections[i].getFromNodeNumber()+1<<" - "<<"Node #"<<connection::connections[i].getToNodeNumber()+1<<" through connection #"<<i+1<<"  =  "<<sum.real()<<" + "<<sum.imagenary()<<"j"<<" OR "<<sum.l<<"ang("<<sum.a<<")"<<endl;
	}
	cout<<endl<<endl;
}

int main()
{
	matrixSolver *mat;
	//getData();
	getDataAwesome();
	mat = evaluateData();
	displayResult(mat);
	cin.ignore();
	cin.get();
	return 0;
}