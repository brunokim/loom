(ns ^{:doc "Defines protocols for graphs, digraphs, and weighted graphs.

Also provides record implementations and constructors for simple graphs --
weighted, unweighted, directed, and undirected. The implementations are based
on adjacency lists."
      :author "Justin Kramer"}
  loom.graph
  (:use [loom.alg-generic :only [bf-traverse]]))

;(set! *warn-on-reflection* true)

;;;
;;; Protocols
;;;

(defprotocol Graph
  (add-nodes* [g nodes] "Add nodes to graph g. See add-nodes")
  (add-edges* [g edges] "Add edges to graph g. See add-edges")
  (remove-nodes* [g nodes] "Remove nodes from graph g. See remove-nodes")
  (remove-edges* [g edges] "Removes edges from graph g. See remove-edges")
  (remove-all [g] "Removes all nodes and edges from graph g")
  (nodes [g] "Return a collection of the nodes in graph g")
  (edges [g] "Edges in g. May return each edge twice in an undirected graph")
  (has-node? [g node] "Return true when node is in g")
  (has-edge? [g n1 n2] "Return true when edge [n1 n2] is in g")
  (neighbors [g] [g node] "Return adjacent nodes, or (partial neighbors g)")
  (degree [g node] "Return the number of nodes adjacent to node"))

(defprotocol Digraph
  (incoming [g node] "Return direct predecessors of node")
  (in-degree [g node] "Return the number direct predecessors to node")
  (transpose [g] "Return a graph with all edges reversed"))

(defprotocol WeightedGraph
  (weight [g] [g n1 n2] "Return weight of edge [n1 n2] or (partial weight g)"))

;; Variadic wrappers

(defn add-nodes
  "Add nodes to graph g. Nodes can be any type of object"
  [g & nodes]
  (add-nodes* g nodes))

(defn add-edges
  "Add edges to graph g. For unweighted graphs, edges take the form [n1 n2].
  For weighted graphs, edges take the form [n1 n2 weight] or [n1 n2], the
  latter defaulting to a weight of 1"
  [g & edges]
  (add-edges* g edges))

(defn remove-nodes
  "Remove nodes from graph g"
  [g & nodes]
  (remove-nodes* g nodes))

(defn remove-edges
  "Remove edges from graph g. Do not include weights"
  [g & edges]
  (remove-edges* g edges))

;;;
;;; Records for simple graphs -- one edge per vertex pair/direction,
;;; loops allowed
;;;
;; TODO: allow custom weight fn?

(defrecord SimpleGraph [nodeset adj])
(defrecord SimpleDigraph [nodeset adj in])
(defrecord SimpleWeightedGraph [nodeset adj])
(defrecord SimpleWeightedDigraph [nodeset adj in])

(def ^{:doc "Weight used when none is given for edges in weighted graphs"}
  *default-weight* 1)

(def default-graph-impls
  {:all
   {:nodes (fn [g]
             (:nodeset g))
    :edges (fn [g]
             (for [n1 (nodes g)
                   n2 (neighbors g n1)]
               [n1 n2]))
    :has-node? (fn [g node]
                 (contains? (:nodeset g) node))
    :has-edge? (fn [g n1 n2]
                 (contains? (get-in g [:adj n1]) n2))
    :degree (fn [g node]
              (count (get-in g [:adj node])))}

   ;; Unweighted graphs store adjacencies as {node #{neighbor}}
   :unweighted
   {:add-nodes* (fn [g nodes]
                  (reduce
                   (fn [g n]
                     (-> g
                         (update-in [:nodeset] conj n)
                         (assoc-in [:adj n] (or ((:adj g) n) #{}))))
                   g nodes))
    :neighbors (fn
                 ([g] (partial neighbors g))
                 ([g node] (get-in g [:adj node])))}
   
   ;; Weighted graphs store adjacencies as {node {neighbor weight}}
   :weighted
   {:add-nodes* (fn [g nodes]
                  (reduce
                   (fn [g n]
                     (-> g
                         (update-in [:nodeset] conj n)
                         (assoc-in [:adj n] (or ((:adj g) n) {}))))
                   g nodes))
    :neighbors (fn
                 ([g] (partial neighbors g))
                 ([g node] (keys (get-in g [:adj node]))))}})

(def default-digraph-impl
  {:incoming (fn [g node]
               (get-in g [:in node]))
   :in-degree (fn [g node]
                (count (get-in g [:in node])))})

(def default-weighted-graph-impl
  {:weight (fn
             ([g] (partial weight g))
             ([g n1 n2] (get-in g [:adj n1 n2])))})

(defn- remove-adj-nodes [m nodes adjacents remove-fn]
  (reduce
   (fn [m n]
     (if (m n)
       (update-in m [n] #(apply remove-fn % nodes))
       m))
   (apply dissoc m nodes)
   adjacents))

(extend SimpleGraph
  Graph
  (assoc (apply merge (map default-graph-impls [:all :unweighted]))

    :add-edges*
    (fn [g edges]
      (reduce
       (fn [g [n1 n2]]
         (-> g
             (update-in [:nodeset] conj n1 n2)
             (update-in [:adj n1] (fnil conj #{}) n2)
             (update-in [:adj n2] (fnil conj #{}) n1)))
       g edges))
    
    :remove-nodes*
    (fn [g nodes]
      (let [nbrs (mapcat #(neighbors g %) nodes)]
        (-> g
            (update-in [:nodeset] #(apply disj % nodes))
            (assoc :adj (remove-adj-nodes (:adj g) nodes nbrs disj)))))
    
    :remove-edges*
    (fn [g edges]
      (reduce
       (fn [g [n1 n2]]
         (-> g
             (update-in [:adj n1] disj n2)
             (update-in [:adj n2] disj n1)))
       g edges))

    :remove-all
    (fn [g]
      (assoc g :nodeset #{} :adj {}))))

(extend SimpleDigraph
  Graph
  (assoc (apply merge (map default-graph-impls [:all :unweighted]))
    
    :add-edges*
    (fn [g edges]
      (reduce
       (fn [g [n1 n2]]
         (-> g
             (update-in [:nodeset] conj n1 n2)
             (update-in [:adj n1] (fnil conj #{}) n2)
             (update-in [:in n2] (fnil conj #{}) n1)))
       g edges))
    
    :remove-nodes*
    (fn [g nodes]
      (let [ins (mapcat #(incoming g %) nodes)
            outs (mapcat #(neighbors g %) nodes)]
        (-> g
            (update-in [:nodeset] #(apply disj % nodes))
            (assoc :adj (remove-adj-nodes (:adj g) nodes ins disj))
            (assoc :in (remove-adj-nodes (:in g) nodes outs disj)))))
    
    :remove-edges*
    (fn [g edges]
      (reduce
       (fn [g [n1 n2]]
         (-> g
             (update-in [:adj n1] disj n2)
             (update-in [:in n2] disj n1)))
       g edges))

    :remove-all
    (fn [g]
      (assoc g :nodeset #{} :adj {} :in {})))

  Digraph
  (assoc default-digraph-impl
    :transpose (fn [g]
                 (assoc g :adj (:in g) :in (:adj g)))))

(extend SimpleWeightedGraph
  Graph
  (assoc (apply merge (map default-graph-impls [:all :weighted]))

    :add-edges*
    (fn [g edges]
      (reduce
       (fn [g [n1 n2 & [w]]]
         (-> g
             (update-in [:nodeset] conj n1 n2)
             (assoc-in [:adj n1 n2] (or w *default-weight*))
             (assoc-in [:adj n2 n1] (or w *default-weight*))))
       g edges))
    
    :remove-nodes*
    (fn [g nodes]
      (let [nbrs (mapcat #(neighbors g %) nodes)]
        (-> g
            (update-in [:nodeset] #(apply disj % nodes))
            (assoc :adj (remove-adj-nodes (:adj g) nodes nbrs dissoc)))))
    
    :remove-edges*
    (fn [g edges]
      (reduce
       (fn [g [n1 n2]]
         (-> g
             (update-in [:adj n1] dissoc n2)
             (update-in [:adj n2] dissoc n1)))
       g edges))

    :remove-all
    (fn [g]
      (assoc g :nodeset #{} :adj {})))
  
  WeightedGraph
  default-weighted-graph-impl)

(extend SimpleWeightedDigraph
  Graph
  (assoc (apply merge (map default-graph-impls [:all :weighted]))
    
    :add-edges*
    (fn [g edges]
      (reduce
       (fn [g [n1 n2 & [w]]]
         (-> g
             (update-in [:nodeset] conj n1 n2)
             (assoc-in [:adj n1 n2] (or w *default-weight*))
             (update-in [:in n2] (fnil conj #{}) n1)))
       g edges))
    
    :remove-nodes*
    (fn [g nodes]
      (let [ins (mapcat #(incoming g %) nodes)
            outs (mapcat #(neighbors g %) nodes)]
        (-> g
            (update-in [:nodeset] #(apply disj % nodes))
            (assoc :adj (remove-adj-nodes (:adj g) nodes ins dissoc))
            (assoc :in (remove-adj-nodes (:in g) nodes outs disj)))))
    
    :remove-edges*
    (fn [g edges]
      (reduce
       (fn [g [n1 n2]]
         (-> g
             (update-in [:adj n1] dissoc n2)
             (update-in [:in n2] disj n1)))
       g edges))

    :remove-all
    (fn [g]
      (assoc g :nodeset #{} :adj {}) :in {}))

  Digraph
  (assoc default-digraph-impl
    :transpose (fn [g]
                 (reduce (fn [tg [n1 n2]]
                             (add-edges* tg [[n2 n1 (weight g n1 n2)]]))
                         (assoc g :adj {} :in {})
                         (edges g))))
  
  WeightedGraph
  default-weighted-graph-impl)

;;;
;;; FlyGraph -- a read-only, ad-hoc graph which uses provided functions to
;;; return values for nodes, edges, etc. Members which are not functions get
;;; returned as-is. Edges can be inferred if nodes and neighbors are provided.
;;; Nodes and edges can be inferred if neighbors and start are provided.
;;;

(defn- call-or-return [f & args]
  (if (or (fn? f)
          (and (instance? clojure.lang.IFn f) (seq args)))
    (apply f args)
    f))

(def ^{:private true} default-flygraph-graph-impl
  {:nodes (fn [g]
            (if (or (:fnodes g) (not (:start g)))
              (call-or-return (:fnodes g))
              (bf-traverse (neighbors g) (:start g))))
   :edges (fn [g]
            (if (:fedges g)
              (call-or-return (:fedges g))
              (for [n (nodes g)
                    nbr (neighbors g n)]
                [n nbr])))
   :neighbors (fn
                ([g] (partial neighbors g))
                ([g node] (call-or-return (:fneighbors g) node)))
   :degree (fn [g node]
             (count (neighbors g node)))})

(def ^{:private true} default-flygraph-digraph-impl
  {:incoming (fn [g node] (call-or-return (:fincoming g) node))
   :in-degree (fn [g node] (count (incoming g node)))})

(def ^{:private true} default-flygraph-weighted-impl
  {:weight (fn [g n1 n2] (call-or-return (:fweight g) n1 n2))})

(defrecord FlyGraph [fnodes fedges fneighbors start])
(defrecord FlyDigraph [fnodes fedges fneighbors fincoming start])
(defrecord WeightedFlyGraph [fnodes fedges fneighbors fweight start])
(defrecord WeightedFlyDigraph [fnodes fedges fneighbors fincoming fweight start])

(extend FlyGraph
  Graph default-flygraph-graph-impl)

(extend FlyDigraph
  Graph default-flygraph-graph-impl
  Digraph default-flygraph-digraph-impl)

(extend WeightedFlyGraph
  Graph default-flygraph-graph-impl
  WeightedGraph default-flygraph-weighted-impl)

(extend WeightedFlyDigraph
  Graph default-flygraph-graph-impl
  Digraph default-flygraph-digraph-impl
  WeightedGraph default-flygraph-weighted-impl)

;;;
;;; Utility functions and constructors
;;;

;; TODO: make this work with read-only graphs?
;; Could also gain speed being impl-specific
(defn subgraph
  "Return a graph without all but the given nodes"
  [g ns]
  (remove-nodes* g (filter (complement (set ns)) (nodes g))))

(defn add-path
  "Add a path of edges connecting the given nodes in order"
  [g & nodes]
  (add-edges* g (partition 2 1 nodes)))

(defn add-cycle
  "Add a cycle of edges connecting the given nodes in order"
  [g & nodes]
  (add-edges* g (partition 2 1 (concat nodes [(first nodes)]))))

(defn graph?
  "Return true if g satisfies the Graph protocol"
  [g]
  (satisfies? Graph g))

(defn directed?
  "Return true if g satisfies the Digraph protocol"
  [g]
  (satisfies? Digraph g))

(defn weighted?
  "Return true if g satisfies the WeightedGraph protocol"
  [g]
  (satisfies? WeightedGraph g))

(defn build-graph
  "Builds up a graph (i.e. adds edges and nodes) from any combination of
  other graphs, adjacency maps, edges, or nodes."
  [g & inits]
  (letfn [(build [g init]
            (cond
             ;; graph
             (graph? init)
             (if (and (weighted? g) (weighted? init))
               (reduce add-edges
                       (add-nodes* g (nodes init))
                       (for [[n1 n2] (edges init)]
                         [n1 n2 (weight init n1 n2)]))
               (-> g
                   (add-nodes* (nodes init))
                   (add-edges* (edges init))))
             ;; adacency map
             (map? init)
             (let [es (if (map? (val (first init)))
                        (for [[n nbrs] init
                              [nbr wt] nbrs]
                          [n nbr wt])
                        (for [[n nbrs] init
                              nbr nbrs]
                          [n nbr]))]
               (-> g
                   (add-nodes* (keys init))
                   (add-edges* es)))
             ;; edge
             (sequential? init) (add-edges g init)
             ;; node
             :else (add-nodes g init)))]
    (reduce build g inits)))

(defn graph
  "Create an unweighted, undirected graph. inits can be edges, adjacency maps,
  or graphs"
  [& inits]
  (apply build-graph (SimpleGraph. #{} {}) inits))

(defn digraph
  "Create an unweighted, directed graph. inits can be edges, adjacency maps,
  or graphs"
  [& inits]
  (apply build-graph (SimpleDigraph. #{} {} {}) inits))

(defn weighted-graph
  [& inits]
  "Create an weighted, undirected graph. inits can be edges, adjacency maps,
  or graphs"
  (apply build-graph (SimpleWeightedGraph. #{} {}) inits))

(defn weighted-digraph
  "Create an weighted, directed graph. inits can be edges, adjacency maps,
  or graphs"
  [& inits]
  (apply build-graph (SimpleWeightedDigraph. #{} {} {}) inits))

(defn fly-graph
  "Create a read-only, ad-hoc graph which uses the provided functions
  to return values for nodes, edges, etc. If any members are not functions,
  they will be returned as-is. Edges can be inferred if nodes and
  neighbors are provided. Nodes and edges can be inferred if neighbors and
  start are provided."
  [& {:keys [nodes edges neighbors incoming weight start]}]
  (cond
   (and incoming weight)
   (WeightedFlyDigraph. nodes edges neighbors incoming weight start)
   incoming
   (FlyDigraph. nodes edges neighbors incoming start)
   weight
   (WeightedFlyGraph. nodes edges neighbors weight start)
   :else
   (FlyGraph. nodes edges neighbors start)))