FROM docker://benjaminkirk/ncar-derecho-wrf:latest
SHELL ["/bin/bash", "-c"]

# Compile WRF
RUN source /container/config_env.sh && \
export WRF_CHEM=0 && \
rm -rf /container/WRF && \
git clone https://github.com/wrf-model/WRF.git /container/WRF && \
cd /container/WRF && \
./clean -a && \
echo "35\n1" | ./configure && \
./compile -j 8 em_real && \
test -f ./main/wrf.exe && \
test -f ./main/real.exe && \
# change permission for all uses to access 
chmod -R a+rx /container

CMD ["/bin/bash"]
