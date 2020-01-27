import numpy as np
import libcdftree
import pytest

def test_insert_sample_no_check():
    libcdftree.init_memory()

    a = np.float32(np.random.uniform(size=(1000,)))
    libcdftree.insert_sample(0, a)

    libcdftree.free_memory()

def test_insert_sample_bad_input_type():
    libcdftree.init_memory()

    a = np.float64(np.random.uniform(size=(1000,)))
    with pytest.raises(RuntimeError):
        libcdftree.insert_sample(0, a)

    libcdftree.free_memory()

def test_insert_sample_bad_input_shape():
    libcdftree.init_memory()

    a = np.float32(np.random.uniform(size=(1000,1)))
    with pytest.raises(RuntimeError):
        libcdftree.insert_sample(0, a)

    libcdftree.free_memory()


def test_samples_to_cdf():
    libcdftree.init_memory()

    data_orig = np.linspace(-4,5,num=20,dtype=np.float32)
    data_copy = np.copy(data_orig)

    libcdftree.insert_sample(0, data_orig)
    v1 = libcdftree.sample_to_cdf(0, data_copy, False)
    v2 = libcdftree.sample_to_cdf(0, data_copy, True)

    assert data_copy is v2
    assert data_copy is not v1
    for a,b in zip(v1, v2):
        assert a == pytest.approx(b)
    for a,b in zip(v1, np.linspace(start=1/20., stop=1, num=20,dtype=np.float32)):
        assert a == pytest.approx(b)

    libcdftree.free_memory()

def test_search_by_cdf():
    libcdftree.init_memory()

    linsp_orig = np.linspace(1/20. - 1/40., 1, num=20, dtype=np.float32)
    linsp_copy = np.linspace(1/20. - 1/40., 1, num=20, dtype=np.float32)

    data_orig = np.linspace(-4,5,num=20,dtype=np.float32)

    libcdftree.insert_sample(0, data_orig)

    v1 = libcdftree.search_element_by_cdf(0, linsp_orig, False)
    v2 = libcdftree.search_element_by_cdf(0, linsp_orig)

    assert v2 is linsp_orig
    assert v1 is not linsp_orig
    for a,b in zip(data_orig, v1):
        print(a,b)
        assert a == pytest.approx(b)
    for a,b in zip(data_orig, v2):
        assert a == pytest.approx(b)

    libcdftree.free_memory()





