#ifndef _CONTACT_FUSION_HPP_
#define _CONTACT_FUSION_HPP_ 1

#include <map>
#include <pair>
#include "logging.hpp"


namespace afidd
{
namespace smv
{

template<typename Subgraph, typename SubgraphKind,
    typename ContactGraph, typename CGVertexToId, typename CGVertexToKind,
    typename TransitionMaker>
class ContactFusion
{
public:

  typedef std::pair<SubgraphKind, typename Subgraph::UserTransitionKey>
      RuleTransition;
  typedef std::map<std::pair<RuleTransition,RuleTransition>,
      TransitionMaker> RuleSet;

  typedef std::pair<SubgraphKind, typename Subgraph::vertex_descriptor>
      VertexRuleTransition;
  typedef std::map<std::pair<VertexRuleTransition,VertexRuleTransition>,
      TransitionMaker> VertexRuleSet;

  Fuse(const std::map<SubgraphKind, Subgraph*>& subgraphs,
    const ContactGraph& cg, const CGVertexToId& cg_vertex_to_id,
    const CGVertexToKind& cg_vertex_to_kind,
    const RuleSet& rule_set) {
    BOOST_LOG_TRIVIAL(debug) << "ContactFusion::Fuse()";



    BOOST_LOG_TRIVIAL(debug) << "ConstactFusion::~Fuse()";
  }

private:
  VertexRuleSet SpecificRules(const RuleSet& rule_set,
    const std::map<SubgraphKind, Subgraph*>& subgraphs) {
    VertexRuleSet vertex_rule_set;
    for (auto& each_rule : rule_set) {
      auto& left_match=each_rule.first.first;
      auto& right_match=each_rule.first.second;
      auto& transition_maker=each_rule.second;
      
    }
  }
};
	
}
}

#endif
