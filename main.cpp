#include <iostream>
#include <cmath>
#include <exception>
using namespace std;

struct polarNum;
class connection;
class node;

polarNum polarParse(float,float);

struct polarNum 
{
	float l;
	float a;
	float real()
	{
		return l*cos(a);
	}
	float imagenary()
	{
		return l*sin(a);
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
	node();
	~node();
	void setConnectionCount(int,int);
	void addConnectionDetails(connection*);
	connection* getNotTraversedConnection();
	connection* getConnectionTo(node*);
	connection* getNextConnection();
	int getSelfNum();
	void resetConnectionCounter();
};

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
	connection();
	~connection();
	void setDetails(node*,node* ,polarNum,polarNum,polarNum);
	void printData();
	void setTraversed(node*,int);
	bool isTraversed();
	void setCurrent(polarNum[],int);
	node* getOtherNode(node*);
};

struct currentPath
{
	node *through;
	connection *from;
	connection *to;
	currentPath *next;
	currentPath *supposeNext;
	int selfNum;
	static int pathCountActual;
	static int maxLoops;
	static int pathCountSuppose;
	static int currentCounter;
	currentPath(node*,currentPath*);
	static currentPath* make(node*,currentPath*);
	static currentPath* start;
	static currentPath* startNewTraversePath(node*);
	bool getNextNode();
	void solidifyPath();
	void convertSupposeToActual();
	void getEquation(polarNum[]);
	bool hasOccured(node*);
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
					cout<<"Infinite Solutions are possible";
					break;
				}
				if(j==variableCount && matA[i][variableCount]!=0)
				{
					cout<<"No solution is possible";
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

connection* node::getConnectionTo(node* _node)
{
	for (int i = 0; i < connectionCount; ++i)
	{
		if(connectionArray[i]->getOtherNode(this) == _node)
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

void connection::setCurrent(polarNum currentArray[],int currentNumber)
{
	int direct;
	if(variable[currentNumber]==-1)
		direct = -1;
	else
		direct = 1;
	for (int i = 0; i < currentPath::maxLoops; ++i)
		currentArray[i]=currentArray[i]+impedence*variable[i]*float(direct);
	currentArray[currentPath::maxLoops] = currentArray[currentPath::maxLoops]+voltage*variable[currentNumber];
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
	if(_prev == NULL)
	{
		start = this;
		pathCountActual = 0;
		pathCountSuppose = 1;
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
		cout<<temp->through->getSelfNum()+1<<" > ";
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
	printActualCurrentPath();
	currentPath* temp = start;
	while(temp->next!=NULL)
	{
		temp->to = temp->through->getConnectionTo(temp->next->through);
		temp->to->setTraversed(temp->through,selfNum);
		temp = temp->next;
	}
	temp->to = temp->through->getConnectionTo(start->through);
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
	if(pathCountSuppose > 2 && (temp=this->through->getConnectionTo(start->through))!=NULL)
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

void currentPath::getEquation(polarNum arr[])
{
	currentPath *temp = this;
	while(temp->next!=this)
	{
		temp->to->setCurrent(arr,this->selfNum);
		temp = temp->next;
	}
	temp->to->setCurrent(arr,this->selfNum);
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


int main()
{	
	node *nodes;
	connection *connections;
	currentPath **currents,*one,*two;
	int nodeCount=0,connectionCount=0,temp,nLoops;
	
	cout<<"Enter the number of nodes : ";
	cin>>nodeCount;
	nodes = new node[nodeCount];
	
	cout<<"Enter the number of connections per node\n";
	for(int i=0;i<nodeCount;i++)	//setting node deatils
	{
		cout<<"Node #"<<i+1<<" : ";
		cin>>temp;
		nodes[i].setConnectionCount(temp,i);
		connectionCount+=temp;
	}
	connectionCount/=2;
	connections = new connection[connectionCount];

	cout<<"Enter the connection details as 'nodeA nodeB Resistance Inductance Capacitance VoltageSourcePeak VoltageSourcePhase CurrentSourcePeak CurrentSourcePhase'\n";
	for(int i=0;i<connectionCount;i++)	//setting connections details
	{
		int nodeA,nodeB;
		float Resistance,Inductance,Capacitance,VoltageSourcePeak,VoltageSourcePhase,CurrentSourcePeak,CurrentSourcePhase;
		cout<<"Connection #"<<i+1<<" : ";
		cin>>nodeA>>nodeB>>Resistance>>Inductance>>Capacitance>>VoltageSourcePeak>>VoltageSourcePhase>>CurrentSourcePeak>>CurrentSourcePhase;
		connections[i].setDetails(&nodes[nodeA-1],&nodes[nodeB-1],polarParse(Resistance,Inductance-Capacitance),polarNum({VoltageSourcePeak,VoltageSourcePhase/180*M_PI}),polarNum({CurrentSourcePeak,CurrentSourcePhase/180*M_PI}));
		nodes[nodeA-1].addConnectionDetails(&connections[i]);
		nodes[nodeB-1].addConnectionDetails(&connections[i]);
	}
	nLoops = connectionCount - nodeCount + 1;
	currents = new currentPath*[nLoops];
	currentPath::maxLoops = nLoops;
	for (int i = 0; i < nodeCount; ++i)
	{
		while((one = currentPath::startNewTraversePath(&nodes[i]))!=NULL)
		{
			currents[currentPath::currentCounter++] = one;
			cout<<"final path:";
			two = one;
			do
			{
				cout<<two->through->getSelfNum()+1<<" > ";
				two = two->next;
			} while (two!=one);
			cout<<two->through->getSelfNum()+1<<endl<<endl;
		}
	}

	polarNum arr[nLoops+1];
	matrixSolver mat(nLoops);
	for (int j = 0; j < nLoops; ++j)
	{
		for (int i = 0; i < nLoops+1; ++i)
		{
			arr[i].l=0;
			arr[i].a=0;
		}	
		currents[j]->getEquation(arr);
		for (int i = 0; i < nLoops+1; ++i)	cout<<arr[i].real()<<" + ";
			cout<<endl;
		mat.addEquations(arr);
	}
	mat.solve();
	for (int i = 0; i < nLoops; ++i)
	{
		cout<<"solution current "<<i+1<<" = "<<(mat.getSolutionFor(i)).real()<<" + "<<(mat.getSolutionFor(i)).imagenary()<<"i"<<endl;
	}
	
	system("pause");
	system("cls");
	return 0;
}