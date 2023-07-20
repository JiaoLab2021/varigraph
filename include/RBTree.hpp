#ifndef RBTree_hpp
#define RBTree_hpp
#include <iostream>

using namespace std;

// 枚举定义结点的颜色
enum Colour
{
	RED,
	BLACK
};

// 红黑树结点的定义
template<class T>
struct RBTreeNode
{
	// 三叉链
	RBTreeNode<T>* _left;
	RBTreeNode<T>* _right;
	RBTreeNode<T>* _parent;

	// 存储的数据
	T _data;

	// 结点的颜色
	int _col; // 红/黑

	// 构造函数
	RBTreeNode(const T& data)
		:_left(nullptr)
		, _right(nullptr)
		, _parent(nullptr)
		, _data(data)
		, _col(RED)
	{}
};

// 正向迭代器
template<class T, class Ref, class Ptr>
struct __TreeIterator
{
	typedef RBTreeNode<T> Node; // 结点的类型
	typedef __TreeIterator<T, Ref, Ptr> Self; // 正向迭代器的类型

    // 构造函数
    __TreeIterator(Node* node)
        :_node(node) // 根据所给结点指针构造一个正向迭代器
    {}

	Ref operator*()
	{
		return _node->_data; // 返回结点数据的引用
	}

	Ptr operator->()
	{
		return &_node->_data; // 返回结点数据的指针
	}

	// 判断两个正向迭代器是否不同
	bool operator!=(const Self& s) const
	{
		return _node != s._node; // 判断两个正向迭代器所封装的结点是否是同一个
	}
	// 判断两个正向迭代器是否相同
	bool operator==(const Self& s) const
	{
		return _node == s._node; // 判断两个正向迭代器所封装的结点是否是同一个
	}

	// 前置++
	Self operator++()
	{
		if (_node->_right) // 结点的右子树不为空
		{
			// 寻找该结点右子树当中的最左结点
			Node* left = _node->_right;
			while (left->_left)
			{
				left = left->_left;
			}
			_node = left; // ++后变为该结点
		}
		else // 结点的右子树为空
		{
			// 寻找孩子不在父亲右的祖先
			Node* cur = _node;
			Node* parent = cur->_parent;
			while (parent&&cur == parent->_right)
			{
				cur = parent;
				parent = parent->_parent;
			}
			_node = parent; // ++后变为该结点
		}
		return *this;
	}

	// 前置--
	Self operator--()
	{
		if (_node->_left) // 结点的左子树不为空
		{
			// 寻找该结点左子树当中的最右结点
			Node* right = _node->_left;
			while (right->_right)
			{
				right = right->_right;
			}
			_node = right; // --后变为该结点
		}
		else // 结点的左子树为空
		{
			// 寻找孩子不在父亲左的祖先
			Node* cur = _node;
			Node* parent = cur->_parent;
			while (parent&&cur == parent->_left)
			{
				cur = parent;
				parent = parent->_parent;
			}
			_node = parent; // --后变为该结点
		}
		return *this;
	}

	Node* _node; // 正向迭代器所封装结点的指针
};

// 反向迭代器---迭代器适配器
template<class Iterator>
struct ReverseIterator
{
	typedef ReverseIterator<Iterator> Self; // 反向迭代器的类型
	typedef typename Iterator::reference Ref; // 结点指针的引用
	typedef typename Iterator::pointer Ptr; // 结点指针

	Iterator _it; // 反向迭代器所封装的正向迭代器

	// 构造函数
	ReverseIterator(Iterator it)
		:_it(it) // 根据所给正向迭代器构造一个反向迭代器
	{}

	Ref operator*()
	{
		return *_it; // 通过调用正向迭代器的operator*返回结点数据的引用
	}
	Ptr operator->()
	{
		return _it.operator->(); // 通过调用正向迭代器的operator->返回结点数据的指针
	}

	// 前置++
	Self& operator++()
	{
		--_it; // 调用正向迭代器的前置--
		return *this;
	}
	// 前置--
	Self& operator--()
	{
		++_it; // 调用正向迭代器的前置++
		return *this;
	}

	bool operator!=(const Self& s) const
	{
		return _it != s._it; // 调用正向迭代器的operator!=
	}
	bool operator==(const Self& s) const
	{
		return _it == s._it; // 调用正向迭代器的operator==
	}
};

// 红黑树的实现
template<class K, class T, class KeyOfT>
class RBTree
{
	typedef RBTreeNode<T> Node; // 结点的类型
public:
	typedef __TreeIterator<T, T&, T*> iterator; // 正向迭代器
	typedef ReverseIterator<iterator> reverse_iterator; // 反向迭代器
    size_t nodeNum = 0;

    // 记录图中节点数
    size_t size()
    {
        return nodeNum;
    }

	reverse_iterator rbegin()
	{
		// 寻找最右结点
		Node* right = _root;
		while (right&&right->_right)
		{
			right = right->_right;
		}
		// 返回最右结点的反向迭代器
		return reverse_iterator(iterator(right));
	}
	reverse_iterator rend()
	{
		// 返回由nullptr构造得到的反向迭代器（不严谨）
		return reverse_iterator(iterator(nullptr));
	}

	iterator begin()
	{
		// 寻找最左结点
		Node* left = _root;
		while (left&&left->_left)
		{
			left = left->_left;
		}
		// 返回最左结点的正向迭代器
		return iterator(left);
	}
	iterator end()
	{
		// 返回由nullptr构造得到的正向迭代器（不严谨）
		return iterator(nullptr);
	}
	// 构造函数
	RBTree()
		:_root(nullptr)
	{}

	// 拷贝构造
	RBTree(const RBTree<K, T, KeyOfT>& t)
	{
		_root = _Copy(t._root, nullptr);
	}

	// 赋值运算符重载（现代写法）
	RBTree<K, T, KeyOfT>& operator=(RBTree<K, T, KeyOfT> t)
	{
		swap(_root, t._root);
		return *this; // 支持连续赋值
	}

	// 析构函数
	~RBTree()
	{
		_Destroy(_root);
		_root = nullptr;
	}

	// 查找函数
	iterator Find(const K& key)
	{
		KeyOfT kot;
		Node* cur = _root;
		while (cur)
		{
			if (key < kot(cur->_data)) // key值小于该结点的值
			{
				cur = cur->_left; // 在该结点的左子树当中查找
			}
			else if (key > kot(cur->_data)) // key值大于该结点的值
			{
				cur = cur->_right; // 在该结点的右子树当中查找
			}
			else // 找到了目标结点
			{
				return iterator(cur); // 返回该结点
			}
		}
		return end(); // 查找失败
	}

	// 插入函数
	pair<iterator, bool> Insert(const T& data)
	{
		if (_root == nullptr) // 若红黑树为空树，则插入结点直接作为根结点
		{
			_root = new Node(data);
			_root->_col = BLACK; // 根结点必须是黑色
            nodeNum++; // 插入成功，则节点数+1
			return make_pair(iterator(_root), true); // 插入成功
		}
		// 1、按二叉搜索树的插入方法，找到待插入位置
		KeyOfT kot;
		Node* cur = _root;
		Node* parent = nullptr;
		while (cur)
		{
			if (kot(data) < kot(cur->_data)) // 待插入结点的key值小于当前结点的key值
			{
				// 往该结点的左子树走
				parent = cur;
				cur = cur->_left;
			}
			else if (kot(data) > kot(cur->_data)) // 待插入结点的key值大于当前结点的key值
			{
				// 往该结点的右子树走
				parent = cur;
				cur = cur->_right;
			}
			else // 待插入结点的key值等于当前结点的key值
			{
				return make_pair(iterator(cur), false); // 插入失败
			}
		}

		// 2、将待插入结点插入到树中
		cur = new Node(data); // 根据所给值构造一个结点
		Node* newnode = cur; // 记录新插入的结点（便于后序返回）
		if (kot(data) < kot(parent->_data)) // 新结点的key值小于parent的key值
		{
			// 插入到parent的左边
			parent->_left = cur;
			cur->_parent = parent;
		}
		else // 新结点的key值大于parent的key值
		{
			// 插入到parent的右边
			parent->_right = cur;
			cur->_parent = parent;
		}

		// 3、若插入结点的父结点是红色的，则需要对红黑树进行调整
		while (parent&&parent->_col == RED)
		{
			Node* grandfather = parent->_parent; // parent是红色，则其父结点一定存在
			if (parent == grandfather->_left) // parent是grandfather的左孩子
			{
				Node* uncle = grandfather->_right; // uncle是grandfather的右孩子
				if (uncle&&uncle->_col == RED) // 情况1：uncle存在且为红
				{
					// 颜色调整
					parent->_col = uncle->_col = BLACK;
					grandfather->_col = RED;

					// 继续往上处理
					cur = grandfather;
					parent = cur->_parent;
				}
				else // 情况2+情况3：uncle不存在 + uncle存在且为黑
				{
					if (cur == parent->_left)
					{
						RotateR(grandfather); // 右单旋

						// 颜色调整
						grandfather->_col = RED;
						parent->_col = BLACK;
					}
					else // cur == parent->_right
					{
						RotateLR(grandfather); // 左右双旋

						// 颜色调整
						grandfather->_col = RED;
						cur->_col = BLACK;
					}
					break; // 子树旋转后，该子树的根变成了黑色，无需继续往上进行处理
				}
			}
			else // parent是grandfather的右孩子
			{
				Node* uncle = grandfather->_left; // uncle是grandfather的左孩子
				if (uncle&&uncle->_col == RED) // 情况1：uncle存在且为红
				{
					// 颜色调整
					uncle->_col = parent->_col = BLACK;
					grandfather->_col = RED;

					// 继续往上处理
					cur = grandfather;
					parent = cur->_parent;
				}
				else // 情况2+情况3：uncle不存在 + uncle存在且为黑
				{
					if (cur == parent->_left)
					{
						RotateRL(grandfather); // 右左双旋

						// 颜色调整
						cur->_col = BLACK;
						grandfather->_col = RED;
					}
					else // cur == parent->_right
					{
						RotateL(grandfather); // 左单旋

						// 颜色调整
						grandfather->_col = RED;
						parent->_col = BLACK;
					}
					break; // 子树旋转后，该子树的根变成了黑色，无需继续往上进行处理
				}
			}
		}
		_root->_col = BLACK; // 根结点的颜色为黑色（可能被情况一变成了红色，需要变回黑色）
        nodeNum++; // 插入成功，则节点数+1
		return make_pair(iterator(newnode), true); // 插入成功
	}

	// 删除函数
	bool Erase(const K& key)
	{
		KeyOfT kot;
		// 用于遍历二叉树
		Node* parent = nullptr;
		Node* cur = _root;
		// 用于标记实际的待删除结点及其父结点
		Node* delParentPos = nullptr;
		Node* delPos = nullptr;
		while (cur)
		{
			if (key < kot(cur->_data)) // 所给key值小于当前结点的key值
			{
				// 往该结点的左子树走
				parent = cur;
				cur = cur->_left;
			}
			else if (key > kot(cur->_data)) // 所给key值大于当前结点的key值
			{
				// 往该结点的右子树走
				parent = cur;
				cur = cur->_right;
			}
			else // 找到了待删除结点
			{
				if (cur->_left == nullptr) // 待删除结点的左子树为空
				{
					if (cur == _root) // 待删除结点是根结点
					{
						_root = _root->_right; // 让根结点的右子树作为新的根结点
						if (_root)
						{
							_root->_parent = nullptr;
							_root->_col = BLACK; // 根结点为黑色
						}
						delete cur; // 删除原根结点
                        nodeNum--; // 删除成功，则节点数-1
						return true;
					}
					else
					{
						delParentPos = parent; // 标记实际删除结点的父结点
						delPos = cur; // 标记实际删除的结点
					}
					break; // 进行红黑树的调整以及结点的实际删除
				}
				else if (cur->_right == nullptr) // 待删除结点的右子树为空
				{
					if (cur == _root) // 待删除结点是根结点
					{
						_root = _root->_left; // 让根结点的左子树作为新的根结点
						if (_root)
						{
							_root->_parent = nullptr;
							_root->_col = BLACK; // 根结点为黑色
						}
						delete cur; // 删除原根结点
                        nodeNum--; // 删除成功，则节点数-1
						return true;
					}
					else
					{
						delParentPos = parent; // 标记实际删除结点的父结点
						delPos = cur; // 标记实际删除的结点
					}
					break; // 进行红黑树的调整以及结点的实际删除
				}
				else // 待删除结点的左右子树均不为空
				{
					// 替换法删除
					// 寻找待删除结点右子树当中key值最小的结点作为实际删除结点
					Node* minParent = cur;
					Node* minRight = cur->_right;
					while (minRight->_left)
					{
						minParent = minRight;
						minRight = minRight->_left;
					}
					cur->_data = minRight->_data; // 将待删除结点的_data改为minRight的_data
					delParentPos = minParent; // 标记实际删除结点的父结点
					delPos = minRight; // 标记实际删除的结点
					break; // 进行红黑树的调整以及结点的实际删除
				}
			}
		}
		if (delPos == nullptr) // delPos没有被修改过，说明没有找到待删除结点
		{
			return false;
		}

		// 记录待删除结点及其父结点（用于后续实际删除）
		Node* del = delPos;
		Node* delP = delParentPos;

		// 调整红黑树
		if (delPos->_col == BLACK) // 删除的是黑色结点
		{
			if (delPos->_left) // 待删除结点有一个红色的左孩子（不可能是黑色）
			{
				delPos->_left->_col = BLACK; // 将这个红色的左孩子变黑即可
			}
			else if (delPos->_right) // 待删除结点有一个红色的右孩子（不可能是黑色）
			{
				delPos->_right->_col = BLACK; // 将这个红色的右孩子变黑即可
			}
			else // 待删除结点的左右均为空
			{
				while (delPos != _root) // 可能一直调整到根结点
				{
					if (delPos == delParentPos->_left) // 待删除结点是其父结点的左孩子
					{
						Node* brother = delParentPos->_right; // 兄弟结点是其父结点的右孩子
						// 情况一：brother为红色
						if (brother->_col == RED)
						{
							delParentPos->_col = RED;
							brother->_col = BLACK;
							RotateL(delParentPos);
							// 需要继续处理
							brother = delParentPos->_right; // 更新brother（否则在本循环中执行其他情况的代码会出错）
						}
						// 情况二：brother为黑色，且其左右孩子都是黑色结点或为空
						if (((brother->_left == nullptr) || (brother->_left->_col == BLACK))
							&& ((brother->_right == nullptr) || (brother->_right->_col == BLACK)))
						{
							brother->_col = RED;
							if (delParentPos->_col == RED)
							{
								delParentPos->_col = BLACK;
								break;
							}
							// 需要继续处理
							delPos = delParentPos;
							delParentPos = delPos->_parent;
						}
						else
						{
							// 情况三：brother为黑色，且其左孩子是红色结点，右孩子是黑色结点或为空
							if ((brother->_right == nullptr) || (brother->_right->_col == BLACK))
							{
								brother->_left->_col = BLACK;
								brother->_col = RED;
								RotateR(brother);
								// 需要继续处理
								brother = delParentPos->_right; // 更新brother（否则执行下面情况四的代码会出错）
							}
							// 情况四：brother为黑色，且其右孩子是红色结点
							brother->_col = delParentPos->_col;
							delParentPos->_col = BLACK;
							brother->_right->_col = BLACK;
							RotateL(delParentPos);
							break; // 情况四执行完毕后调整一定结束
						}
					}
					else // delPos == delParentPos->_right // 待删除结点是其父结点的左孩子
					{
						Node* brother = delParentPos->_left; // 兄弟结点是其父结点的左孩子
						// 情况一：brother为红色
						if (brother->_col == RED) // brother为红色
						{
							delParentPos->_col = RED;
							brother->_col = BLACK;
							RotateR(delParentPos);
							// 需要继续处理
							brother = delParentPos->_left; // 更新brother（否则在本循环中执行其他情况的代码会出错）
						}
						// 情况二：brother为黑色，且其左右孩子都是黑色结点或为空
						if (((brother->_left == nullptr) || (brother->_left->_col == BLACK))
							&& ((brother->_right == nullptr) || (brother->_right->_col == BLACK)))
						{
							brother->_col = RED;
							if (delParentPos->_col == RED)
							{
								delParentPos->_col = BLACK;
								break;
							}
							// 需要继续处理
							delPos = delParentPos;
							delParentPos = delPos->_parent;
						}
						else
						{
							// 情况三：brother为黑色，且其右孩子是红色结点，左孩子是黑色结点或为空
							if ((brother->_left == nullptr) || (brother->_left->_col == BLACK))
							{
								brother->_right->_col = BLACK;
								brother->_col = RED;
								RotateL(brother);
								// 需要继续处理
								brother = delParentPos->_left; // 更新brother（否则执行下面情况四的代码会出错）
							}
							// 情况四：brother为黑色，且其左孩子是红色结点
							brother->_col = delParentPos->_col;
							delParentPos->_col = BLACK;
							brother->_left->_col = BLACK;
							RotateR(delParentPos);
							break; // 情况四执行完毕后调整一定结束
						}
					}
				}
			}
		}
		// 进行实际删除
		if (del->_left == nullptr) // 实际删除结点的左子树为空
		{
			if (del == delP->_left) // 实际删除结点是其父结点的左孩子
			{
				delP->_left = del->_right;
				if (del->_right)
					del->_right->_parent = delP;
			}
			else // 实际删除结点是其父结点的右孩子
			{
				delP->_right = del->_right;
				if (del->_right)
					del->_right->_parent = delP;
			}
		}
		else // 实际删除结点的右子树为空
		{
			if (del == delP->_left) // 实际删除结点是其父结点的左孩子
			{
				delP->_left = del->_left;
				if (del->_left)
					del->_left->_parent = delP;
			}
			else // 实际删除结点是其父结点的右孩子
			{
				delP->_right = del->_left;
				if (del->_left)
					del->_left->_parent = delP;
			}
		}
		delete del; // 实际删除结点
        nodeNum--; // 删除成功，则节点数-1
		return true;
	}

private:
	// 拷贝树
	Node* _Copy(Node* root, Node* parent)
	{
		if (root == nullptr)
		{
			return nullptr;
		}
		Node* copyNode = new Node(root->_data);
		copyNode->_parent = parent;
		copyNode->_left = _Copy(root->_left, copyNode);
		copyNode->_right = _Copy(root->_right, copyNode);
		return copyNode;
	}

	// 析构函数子函数
	void _Destroy(Node* root)
	{
		if (root == nullptr)
		{
			return;
		}
		_Destroy(root->_left);
		_Destroy(root->_right);
		delete root;
	}

	// 左单旋
	void RotateL(Node* parent)
	{
		Node* subR = parent->_right;
		Node* subRL = subR->_left;
		Node* parentParent = parent->_parent;

		// 建立subRL与parent之间的联系
		parent->_right = subRL;
		if (subRL)
			subRL->_parent = parent;

		// 建立parent与subR之间的联系
		subR->_left = parent;
		parent->_parent = subR;

		// 建立subR与parentParent之间的联系
		if (parentParent == nullptr)
		{
			_root = subR;
			_root->_parent = nullptr;
		}
		else
		{
			if (parent == parentParent->_left)
			{
				parentParent->_left = subR;
			}
			else
			{
				parentParent->_right = subR;
			}
			subR->_parent = parentParent;
		}
	}

	// 右单旋
	void RotateR(Node* parent)
	{
		Node* subL = parent->_left;
		Node* subLR = subL->_right;
		Node* parentParent = parent->_parent;

		// 建立subLR与parent之间的联系
		parent->_left = subLR;
		if (subLR)
			subLR->_parent = parent;

		// 建立parent与subL之间的联系
		subL->_right = parent;
		parent->_parent = subL;

		// 建立subL与parentParent之间的联系
		if (parentParent == nullptr)
		{
			_root = subL;
			_root->_parent = nullptr;
		}
		else
		{
			if (parent == parentParent->_left)
			{
				parentParent->_left = subL;
			}
			else
			{
				parentParent->_right = subL;
			}
			subL->_parent = parentParent;
		}
	}

	// 左右双旋
	void RotateLR(Node* parent)
	{
		RotateL(parent->_left);
		RotateR(parent);
	}

	// 右左双旋
	void RotateRL(Node* parent)
	{
		RotateR(parent->_right);
		RotateL(parent);
	}

	Node* _root; // 红黑树的根结点
};

namespace cl // 防止命名冲突
{
	template<class K, class V>
	class map
	{
		// 仿函数
		struct MapKeyOfT
		{
			const K& operator()(const pair<K, V>& kv) // 返回键值对当中的键值Key
			{
				return kv.first;
			}
		};
	public:
		typedef typename RBTree<K, pair<K, V>, MapKeyOfT>::iterator iterator; // 正向迭代器
		typedef typename RBTree<K, pair<K, V>, MapKeyOfT>::reverse_iterator reverse_iterator; // 反向迭代器
		
        // 返回图中节点数
        size_t size()
        {
            return _t.size();
        }

		iterator begin()
		{
			return _t.begin();
		}
		iterator end()
		{
			return _t.end();
		}
		
		reverse_iterator rbegin()
		{
			return _t.rbegin();
		}
		reverse_iterator rend()
		{
			return _t.rend();
		}

		// 插入函数
		pair<iterator, bool> insert(const pair<const K, V>& kv)
		{
			return _t.Insert(kv);
		}
		// []运算符重载函数
		V& operator[](const K& key)
		{
			auto ret = insert(make_pair(key, V()));
			iterator it = ret.first;
			return it->second;
		}
		// 删除函数
		void erase(const K& key)
		{
			_t.Erase(key);
		}
		// 查找函数
		iterator find(const K& key)
		{
			return _t.Find(key);
		}
	private:
		RBTree<K, pair<K, V>, MapKeyOfT> _t;
	};
}

#endif