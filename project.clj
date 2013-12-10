(defproject brunokim/loom "0.4.2-SNAPSHOT"
  :description "Graph library for Clojure"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.4.0"]]
  :url "https://github.com/brunokim/loom"
  :profiles {:dev 
             {:dependencies [[org.clojure/clojure "1.4.0"]]}}
  :aliases {"release" ["do" "clean," "with-profile" "default" "deploy" "clojars"]}
  :aot [loom.graph])
