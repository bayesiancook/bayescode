#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

/*
====================================================================================================
  ~*~ Pure interfaces ~*~
==================================================================================================*/
class AnnotatedTree {
  public:
    using NodeIndex = int;
    using TagName = std::string;
    using TagValue = std::string;
    using ChildrenList = std::vector<NodeIndex>;

    virtual const ChildrenList& children(NodeIndex) const = 0;
    virtual NodeIndex parent(NodeIndex) const = 0;
    virtual NodeIndex root() const = 0;
    virtual std::size_t nb_nodes() const = 0;
    virtual TagValue tag(NodeIndex, TagName) const = 0;

    virtual ~AnnotatedTree() = default;
};

class TreeParser {
  public:
    virtual const AnnotatedTree& get_tree() const = 0;
    virtual ~TreeParser() = default;
};

/*
====================================================================================================
  ~*~ Implementations ~*~
==================================================================================================*/
class DoubleListAnnotatedTree : public AnnotatedTree {
  public:
    using Node = std::unordered_map<std::string, std::string>;

    std::vector<Node> nodes_;

    // invariant: same length as nodes
    // invariant: only one node has parent -1
    // element i is the index of parent of i in nodes (-1 for root)
    std::vector<int> parent_;

    // invariant: consistent with parent
    // element i is a vector of indices corresponding to the children of node i
    std::vector<std::vector<int>> children_;

    // invariant: node with index root is only node with parent -1
    NodeIndex root_;

  public:
    const ChildrenList& children(NodeIndex node) const final { return children_.at(node); }

    NodeIndex parent(NodeIndex node) const final { return parent_.at(node); }

    NodeIndex root() const final { return root_; }

    std::size_t nb_nodes() const final { return nodes_.size(); }

    TagValue tag(NodeIndex node, TagName tag) const final {
        if (nodes_.at(node).count(tag) != 0) {
            return nodes_.at(node).at(tag);
        } else {
            return "";
        }
    }
};

/*================================================================================================*/
struct NHXException : public std::runtime_error {
    NHXException(std::string s = "") : std::runtime_error(s) {}
};

/*================================================================================================*/
class NHXParser : public TreeParser {
    // list of tokens for lexer
    enum TokenType {
        OpenParenthesis,
        CloseParenthesis,
        Colon,
        Semicolon,
        Comma,
        Equal,
        NHXOpen,
        CommentOpen,
        BracketClose,
        Identifier,
        Invalid
    };
    std::map<TokenType, std::regex> token_regexes{{OpenParenthesis, std::regex("\\(")},
                                                  {CloseParenthesis, std::regex("\\)")},
                                                  {Colon, std::regex(":")},
                                                  {Semicolon, std::regex(";")},
                                                  {Comma, std::regex(",")},
                                                  {Equal, std::regex("=")},
                                                  {NHXOpen, std::regex("\\[&&NHX:")},
                                                  {CommentOpen, std::regex("\\[")},
                                                  {BracketClose, std::regex("\\]")},
                                                  {Identifier, std::regex("[a-zA-Z0-9._-]+")}};
    using Token = std::pair<TokenType, std::string>;  // first: index of token, second: token value

    // input/output
    DoubleListAnnotatedTree tree;

    // state during parsing
    using scit = std::string::const_iterator;
    scit it;
    Token next_token{Invalid, ""};
    std::string input{""};
    int next_node{0};

    [[noreturn]] void error(std::string s) {
        std::stringstream ss;
        ss << s;
        ss << "Error at position " << std::distance(scit(input.begin()), it) << ":\n";
        bool at_begining = it - 15 <= input.begin();
        bool at_end = it + 15 >= input.end();
        ss << "\t" << (at_begining ? "" : "...")
           << std::string(at_begining ? input.begin() : it - 15, at_end ? input.end() : it + 15)
           << (at_end ? "" : "...") << "\n";
        ss << "\t" << (at_begining ? "" : "   ")
           << std::string(at_begining ? std::distance(scit(input.begin()), it) : 15, ' ') << "^\n";
        throw NHXException(ss.str());
    }

    std::string expect(TokenType type) {
        find_token();
        if (next_token.first != type) {
            error("Error: expected token " + std::to_string(type) + " but got token " +
                  std::to_string(next_token.first) + "(" + next_token.second + ") instead.\n");
        } else {
            return next_token.second;
        }
    }

    // lexer
    void find_token() {
        while (std::isspace(*it)) {
            it++;
        }
        int token_number{0};
        for (auto token_regex : token_regexes) {
            std::smatch m;
            if (std::regex_search(it, it + 64, m, token_regex.second) and m.prefix() == "") {
                if (token_regex.first == CommentOpen) {  // support of comments
                    std::string comment_close{"]"};
                    it = std::search(it, scit(input.end()), comment_close.begin(),
                                     comment_close.end()) +
                         1;
                    find_token();
                    return;
                }
                next_token = Token{token_regex.first, m[0]};
                it += std::string(m[0]).size();
                // std::cout << "found token " << m[0] << std::endl;
                return;
            }
            token_number++;
        }
        next_token = Token{Invalid, std::string(it, it + 10) + "..."};  // no token found in chain
    }

    // parser
    void node_nothing(int number, int parent) {
        tree.nodes_.emplace_back();
        tree.parent_.push_back(parent);
        tree.children_.emplace_back();
        if (parent != -1) {
            tree.children_.at(parent).push_back(number);
        }

        find_token();
        switch (next_token.first) {
            case Identifier:
                tree.nodes_[number]["name"] = next_token.second;
                node_name(number, parent);
                break;
            case Colon:
                node_length(number, parent);
                break;
            case NHXOpen:
                data(number, parent);
                break;
            case OpenParenthesis:
                next_node++;
                node_nothing(next_node, number);
                break;
            default:
                node_end(parent);
        }
    }

    void node_name(int number, int parent) {
        find_token();
        switch (next_token.first) {
            case Colon:
                node_length(number, parent);
                break;
            case NHXOpen:
                data(number, parent);
                break;
            default:
                node_end(parent);
        }
    }

    void node_length(int number, int parent) {
        tree.nodes_[number]["length"] = expect(Identifier);

        find_token();
        switch (next_token.first) {
            case NHXOpen:
                data(number, parent);
                break;
            default:
                node_end(parent);
        }
    }

    void node_end(int parent) {
        switch (next_token.first) {
            case Comma:
                next_node++;
                node_nothing(next_node, parent);
                break;
            case CloseParenthesis: {
                if (parent != -1) {
                    node_name(parent, tree.parent_.at(parent));
                }
                break;
            }
            case Semicolon:
                break;
            default:
                error("Error: unexpected token " + next_token.second + '\n');
        }
    }

    void data(int number, int parent) {
        find_token();
        if (next_token.first == BracketClose) {
            find_token();
            node_end(parent);
        } else if (next_token.first == Identifier) {
            std::string tag = next_token.second;
            expect(Equal);
            tree.nodes_[number][tag] = expect(Identifier);
            data(number, parent);
        } else if (next_token.first == Colon) {
            data(number, parent);
        } else {
            error("Error: improperly formatted contents in NHX data. Found unexpected token " +
                  next_token.second + '\n');
        }
    }

  public:
    NHXParser(std::istream& is) {
        tree = DoubleListAnnotatedTree();
        tree.root_ = 0;
        input = std::string(std::istreambuf_iterator<char>(is), {});
        it = input.begin();
        node_nothing(0, -1);
    }

    const AnnotatedTree& get_tree() const final { return tree; }
};
