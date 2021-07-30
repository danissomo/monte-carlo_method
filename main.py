import PyKDL as kdl
import time
import math
import kdl_parser_py.urdf as urdf
import random
import numpy as np


def get_chain(filename, base_link="base_link", end_link="ee_link"):
    ok, tree = urdf.treeFromFile(filename)
    if ok:
        print("Info: Loading ", filename, " completed")
    else:
        print("Error: Loading ", filename, "failed")
        exit(1)
    chain = tree.getChain(base_link, end_link)
    return chain


def gen_rand_buff(size):
    buf = kdl.JntArray(size)
    for i in range(size):
        buf[i] = random.uniform(-3, 3)
    return buf


def calculate_all(dyn_model, size, qd_prev, dt):
    q = gen_rand_buff(size)
    qd = gen_rand_buff(size)
    grav_torques = kdl.JntArray(size)
    mass_matrix = kdl.JntSpaceInertiaMatrix(size)
    coriolis_vector = kdl.JntArray(size)
    dyn_model.JntToGravity(q, grav_torques)
    dyn_model.JntToMass(q, mass_matrix)
    dyn_model.JntToCoriolis(q, qd, coriolis_vector)
    qdd = kdl.JntArray(size)
    for i in range(size):
        qdd[i] = (qd[i] - qd_prev[i]) / dt
    mass_vec = kdl.JntArray(size)
    kdl.Multiply(mass_matrix, qdd, mass_vec)
    
    return q, qd, qdd, grav_torques, mass_vec, coriolis_vector


def csv_vec(name, size):
    buf = ""
    for i in range(size):
        buf+= name + i.__str__()+ ","
    return buf


def write_results_to_file(filename, time_stamp=[],  q=[], qd=[], qdd=[], grav_torques=[], inertia_torques=[], coriolis=[]):
    size = len(time_stamp)
    n = q[0].rows()
    file = open(filename, 'w')
    file.write("i,"+ "time,"+ csv_vec("q", n)+"," + csv_vec("qd", n)+"," + csv_vec("qdd", n)+"," + csv_vec("grav_trq", n)+"," + csv_vec("inert", n)+"," + csv_vec("coriolis", n)+ f"\n")
    for i in range(size):
        file.write(i.__str__() + "," + time_stamp[i].__str__() + ",")
        for j in range(n):
            file.write("{:.18f}".format(q[i][j]) + ",")
        file.write(",")
        for j in range(n):
            file.write("{:.18f}".format(qd[i][j]) + ",")
        file.write(",")
        for j in range(n):
            file.write("{:.18f}".format(qdd[i][j]) + ",")
        file.write(",")
        for j in range(n):
            file.write("{:.18f}".format(grav_torques[i][j]) + ",")
        file.write(",")
        for j in range(n):
            file.write("{:.18f}".format(inertia_torques[i][j]) + ",")
        file.write(",")
        for j in range(n):
            file.write("{:.18f}".format(coriolis[i][j]) + ",")
        file.write(f"\n")
    file.close()
def calc_delta_write(filename, time_stamp=[],  q=[], qd=[], qdd=[], grav_torques=[], inertia_torques=[], coriolis=[]):
    size = len(time_stamp)
    n = q[0].rows()
    file = open(filename, 'w')
    file.write("i,"+ "time,"+ csv_vec("d_q", n)+"," + csv_vec("d_qd", n)+"," + csv_vec("d_qdd", n)+"," + csv_vec("d_grav_trq", n)+"," + csv_vec("d_inert", n)+"," + csv_vec("d_coriolis", n)+ f"\n")
    file.write("0," + time_stamp[0].__str__() + ",")
    for j in range(n):
        file.write("0,")
    file.write(",")
    for j in range(n):
        file.write("0,")
    file.write(",")
    for j in range(n):
        file.write("0,")
    file.write(",")
    for j in range(n):
        file.write("0,")
    file.write(",")
    for j in range(n):
        file.write("0,")
    file.write(",")
    for j in range(n):
        file.write("0,")
    file.write(f"\n")
    
    for i in range(size):
        file.write(i.__str__() + "," + time_stamp[i].__str__() + ",")
        for j in range(n):
            file.write("{:.18f}".format(q[i][j] - q[i-1][j]) + ",")
        file.write(",")
        for j in range(n):
            file.write("{:.18f}".format(qd[i][j]-qd[i-1][j]) + ",")
        file.write(",")
        for j in range(n):
            file.write("{:.18f}".format(qdd[i][j]-qdd[i-1][j]) + ",")
        file.write(",")
        for j in range(n):
            file.write("{:.18f}".format(grav_torques[i][j]-grav_torques[i-1][j]) + ",")
        file.write(",")
        for j in range(n):
            file.write("{:.18f}".format(inertia_torques[i][j]-inertia_torques[i-1][j]) + ",")
        file.write(",")
        for j in range(n):
            file.write("{:.18f}".format(coriolis[i][j]) + ",")
        file.write(f"\n")
def main():
    dt = 0.008
    test_step_count = 10000
    chain = get_chain("ur5.xml")
    n = 6

    grav = kdl.Vector(0, 0, -9.82)
    dyn_model = kdl.ChainDynParam(chain, grav)
    q = []
    qd = []
    qdd = []
    grav_torques = []
    mass_vec = []
    coriolis_vec = []
    time_stamp = []
    start_time = time.time()
    buf_q, buf_qd, buf_qdd, buf_grav_torques, buf_mass_vec, buf_coriolis_vec = calculate_all(
        dyn_model, n, kdl.JntArray(n), dt)
    q.append(buf_q) 
    qd.append(buf_qd)
    qdd.append(buf_qdd)
    grav_torques.append(buf_grav_torques)
    mass_vec.append(buf_mass_vec)
    coriolis_vec.append(buf_coriolis_vec)
    time_stamp.append(time.time())
    
    for i in range(1, test_step_count):
        buf_q, buf_qd, buf_qdd, buf_grav_torques, buf_mass_vec, buf_coriolis_vec = calculate_all(
        dyn_model, n, qd[i-1], dt)
        q.append(buf_q) 
        qd.append(buf_qd)
        qdd.append(buf_qdd)
        grav_torques.append(buf_grav_torques)
        mass_vec.append(buf_mass_vec)
        coriolis_vec.append(buf_coriolis_vec)
        time_stamp.append(time.time()-start_time)
    write_results_to_file("test.csv", time_stamp, q, qd, qdd, grav_torques, mass_vec, coriolis_vec)
    calc_delta_write("delta_test.csv", time_stamp, q, qd, qdd, grav_torques, mass_vec, coriolis_vec)
main()
