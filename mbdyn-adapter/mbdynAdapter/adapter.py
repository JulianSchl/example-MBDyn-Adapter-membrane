import numpy as np
from mpi4py import MPI
import precice

class MBDynAdapter:
    def __init__(self, mbdHelper, configFileName = "../precice-config.xml"):
        self.mbd = mbdHelper
        self.interface = precice.Interface("Structure_Solver",configFileName , 0, 1) # proc no, nprocs
        self.dim = self.interface.get_dimensions()
        nodes = self.mbd.getNodes()

        # if self.dim == 2:
            # self.nodesid = np.where(nodes[:,2] < 1e-6)
            # nodes = nodes[self.nodesid]

        self.nnodes = len(nodes)

        nmeshID = self.interface.get_mesh_id("Structure_Nodes")
        self.nodeVertexIDs = self.nnodes*[0.0]
        #self.interface.set_mesh_vertices(nmeshID, self.nnodes, np.ravel(nodes[:,:self.dim]), self.nodeVertexIDs)
        #self.interface.set_mesh_vertices(nmeshID, np.ravel(nodes[:,:self.dim]))
        self.interface.set_mesh_vertices(nmeshID, np.array(nodes[:,:self.dim]))
        #self.displacements = np.array(self.dim*self.nnodes*[0.0])
        self.displacements = np.zeros((self.nnodes,self.dim))

        # ccs = self.mbd.getCellCenters()
        # self.ncells = len(ccs)
        # cmeshID = self.interface.get_mesh_id("Structure_CellCenters")
        # self.cellVertexIDs = self.ncells*[0.0]
        # self.interface.set_mesh_vertices(cmeshID, self.ncells, np.ravel(ccs[:,:self.dim]), self.cellVertexIDs)
        self.force = np.array(self.dim*self.nnodes*[0.0])

        self.displacementsID = self.interface.get_data_id("DisplacementDelta", nmeshID)
        self.forceID = self.interface.get_data_id("Force", nmeshID)

        self.dt = self.interface.initialize()
        self.mbd.controlDict["timeStep"] = self.dt
        self.mbd.initializeMBDyn()
        #import ipdb; ipdb.set_trace()


        if (self.interface.is_action_required(precice.action_write_initial_data())):
            print("Writing initial data")
            self.interface.write_block_vector_data(self.displacementsID, self.nodeVertexIDs, self.displacements)
            self.interface.mark_action_fulfilled(precice.action_write_initial_data())

        self.interface.initialize_data()
        if (self.interface.is_read_data_available()):
            print("Reading initial data")
            self.force = self.interface.read_block_vector_data(self.forceID, self.nodeVertexIDs)
        print(4)
        #import ipdb; ipdb.set_trace()

    def runPreCICE(self):
        iteration = 0
        previousDisplacements = self.mbd.getDisplacements()
        while (self.interface.is_coupling_ongoing()):
            print("is write")
            if (self.interface.is_action_required(precice.action_write_iteration_checkpoint())):
                self.interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

            
            #
            if self.interface.is_read_data_available():
            	print("reading")
            	self.force = self.interface.read_block_vector_data(self.forceID, self.nodeVertexIDs)
            #
            
            if self.dim == 2:
                f = np.zeros((self.nnodes,3))
                f[:,:self.dim] = np.reshape(self.force,(-1,self.dim))
            else:
                #import ipdb; ipdb.set_trace()
                f = np.reshape(self.force,(-1,3))
            self.mbd.setForces(f)

            if self.mbd.solve(False):
                break
            displacements = self.mbd.getDisplacements()
            #import ipdb; ipdb.set_trace()
            relDisplacements = displacements - previousDisplacements
            print("write")


            #
            if self.interface.is_write_data_required(self.dt):
            	self.interface.write_block_vector_data(self.displacementsID, self.nodeVertexIDs, np.array(relDisplacements))
            #
            
            
            #import ipdb; ipdb.set_trace()
            '''self.interface.write_block_vector_data(self.displacementsID, self.nodeVertexIDs, np.array(relDisplacements))'''
            print("pre error")
            #import ipdb; ipdb.set_trace()
            self.interface.advance(self.dt)
            print("post error")
            print("read")
            #self.interface.read_block_vector_data(self.forceID, self.nodeVertexIDs)

            if (self.interface.is_action_required(precice.action_read_iteration_checkpoint())): # i.e. not yet
                self.interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())
            else:
                previousDisplacements = displacements.copy()
                iteration += 1

                if self.mbd.solve(True):
                    break
                if iteration % self.mbd.controlDict['output frequency']  == 0:
                    self.mbd.writeVTK(iteration)
            print("loop")
        self.mbd.finalize()
