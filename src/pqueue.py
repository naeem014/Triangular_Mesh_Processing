import math
import numpy

class Pqueue:
    def __init__(self, l_t_size):
        self.list_ = []
        self.l_t = numpy.full(l_t_size, -1, dtype=int)

    def empty(self):
        return len(self.list_) == 0
    
    def size(self):
        return len(self.list_)

    def parent(self, index):
        return math.floor((index - 1) / 2)

    def left(self, index):
        return (2 * index) + 1

    def right(self, index):
        return (2 * index) + 2

    def swap(self, index1, index2):
        self.l_t[self.list_[index1][1]] = index2
        self.l_t[self.list_[index2][1]] = index1
        temp = self.list_[index1]
        self.list_[index1] = self.list_[index2]
        self.list_[index2] = temp
        

    def moveDown(self, index):
        left_index = self.left(index)
        right_index = self.right(index)
        final_index = index

        if (left_index <= self.size() - 1 and self.list_[left_index][0] < self.list_[final_index][0]):
            final_index = left_index
        if (right_index <= self.size() - 1 and self.list_[right_index][0] < self.list_[final_index][0]):
            final_index = right_index
        if final_index != index:
            self.swap(final_index, index)
            self.moveDown(final_index)

    def moveUp(self, index):
        while (True):
            parent_index = self.parent(index)
            if parent_index >= 0 and self.list_[parent_index][0] > self.list_[index][0]:
                self.swap(self.parent(index), index)
                index = parent_index
            else:
                break

    def insert(self, item):
        self.list_.append(item)
        self.l_t[item[1]] = self.size() - 1
        self.moveUp(self.size() - 1)
    
    def updatePriority(self, index, p):
        final_index = index
        oldP = self.list_[index][0]
        self.list_[index][0] = p
        if p > oldP:
            self.moveDown(index)
        else:
            self.moveUp(index)

    def pop(self):
        res = self.list_[0]
        self.list_[0] = self.list_[self.size() - 1]
        self.moveDown(0)
        self.list_.pop()
        return res

