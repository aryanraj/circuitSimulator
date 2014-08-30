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
	polarNum operator+(polarNum &b)
	{
		return polarParse(this->l*cos(this->a)+b.l*cos(b.a),this->l*sin(this->a)+b.l*sin(b.a));
	}
	polarNum operator-(polarNum &b)
	{
		return polarParse(this->l*cos(this->a)-b.l*cos(b.a),this->l*sin(this->a)-b.l*sin(b.a));
	}
	polarNum operator/(polarNum &b)
	{
		polarNum temp;
		temp.l = this->l/b.l;
		temp.a = this->a-b.a;
		return temp;
	}
	polarNum operator*(polarNum &b)
	{
		polarNum temp;
		temp.l = this->l*b.l;
		temp.a = this->a+b.a;
		return temp;
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
	polarNum impedence;
	polarNum voltage;
	polarNum current;
public:
	connection();
	~connection();
	void setDetails(node*,node* /*,polarNum,polarNum,polarNum*/);
	void printData();
	void setTraversed();
	bool isTraversed();
	node* getOtherNode(node*);
};

struct currentPath
{
	node *through;
	connection *from;
	connection *to;
	currentPath *next;
	currentPath *supposeNext;
	static int pathCountActual;
	static int pathCountSuppose;
	currentPath(node*,currentPath*);
	static currentPath* make(node*,currentPath*);
	static currentPath* start;
	static currentPath* startNewTraversePath(node*);
	bool getNextNode();
	void solidifyPath();
	void convertSupposeToActual();
	bool hasOccured(node*);
	~currentPath();
	static void printActualCurrentPath();
	static void printSupposeCurrentPath();
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
}

void connection::setDetails(node *prev, node *next/*,polarNum imp = 0,polarNum volt = 0,polarNum cur = 0*/)
{
	this->fromNode = prev;
	this->toNode = next;
	//this->impedence = imp;
	//this->voltage = volt;
	//this->current = cur;
}

void connection::printData()
{
	cout<<impedence.l<<"  "<<impedence.a/M_PI*180<<" : "<<voltage.l<<"  "<<voltage.a/M_PI*180<<" : "<<current.l<<"  "<<current.a/M_PI*180;
}

void connection::setTraversed()
{
	this->traversed=true;
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
	cout<<"current Path : \n";
	while(temp!=NULL)
	{
		cout<<temp->through->getSelfNum()+1<<" > ";
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
		temp->to->setTraversed();
		temp = temp->next;
	}
	temp->to = temp->through->getConnectionTo(start->through);
	temp->to->setTraversed();
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

polarNum polarParse(float a, float b)
{
	polarNum temp;
	temp.l = sqrt(a*a+b*b);
	temp.a = atan(b/a);
	return temp;
}

int main()
{	
	node *nodes;
	connection *connections;
	currentPath *one,*two;
	int nodeCount=0,connectionCount=0,temp;
	
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
		cin>>nodeA>>nodeB/*>>Resistance>>Inductance>>Capacitance>>VoltageSourcePeak>>VoltageSourcePhase>>CurrentSourcePeak>>CurrentSourcePhase*/;
		connections[i].setDetails(&nodes[nodeA-1],&nodes[nodeB-1]/*,polarParse(Resistance,Inductance-Capacitance),polarNum({VoltageSourcePeak,VoltageSourcePhase/180*M_PI}),polarNum({CurrentSourcePeak,CurrentSourcePhase/180*M_PI})*/);
		nodes[nodeA-1].addConnectionDetails(&connections[i]);
		nodes[nodeB-1].addConnectionDetails(&connections[i]);
	}
	for (int i = 0; i < nodeCount; ++i)
	{
		while((one = currentPath::startNewTraversePath(&nodes[i]))!=NULL)
		{
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

	system("pause");
	system("cls");
	return 0;
}