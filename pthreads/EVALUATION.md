# Experimental Evaluation

Tools:

- [https://sbatchman.readthedocs.io/en/latest/](https://sbatchman.readthedocs.io/en/latest/)
- [https://github.com/ThomasPasquali/MtxMan](https://github.com/ThomasPasquali/MtxMan)

Clone the repo using:

```bash
git clone --bare git@github.com:Sasso0101/thesis.git BFS
cd BFS

git worktree add main
cd main
git submodule init
git submodule update

cd ..
git worktree add cnalock
cd cnalock
git submodule init
git submodule update

cd ..
git worktree add mcslock
cd mcslock
git submodule init
git submodule update

cd ..
git worktree add atomic
cd atomic
git submodule init
git submodule update

cd ..
git worktree add openmp
cd openmp
git submodule init
git submodule update
```

1) Init SbatchMan in the root folder (`BFS`)
2) Download datasets with MtxMan into the root folder (`BFS/datasets`)
3) From the root folder (`BFS`), compile with `cp main/scripts/compile_hca.sh . && ./compile_hca.sh` 
3) Each worktree will have its own yaml files. Launch them from their folder (e.g., for main from `BFS/main`)