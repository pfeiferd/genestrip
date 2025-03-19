package org.metagene.genestrip.util;

import java.util.Collection;
import java.util.Iterator;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;

/**
 * A super simple implementation of BlockingQueue - just good enough in the context of genestrip and no more.
 * Still way more efficient than the standard implementation from the JDK for this purpose.
 */
public class SimpleBlockinqQueue<T> implements BlockingQueue<T> {
    private final Object lock = new Object();

    private final T[] array;
    private final int capacity;
    private int size = 0;
    private int head = 0;
    private int tail = 0;

    @SuppressWarnings("unchecked")
    public SimpleBlockinqQueue(int capacity) {
        this.capacity = capacity;
        array = (T[]) new Object[capacity];
    }


    @Override
    public final void put(final T item) throws InterruptedException {
        synchronized (lock) {
            while (size == capacity) {
                lock.wait();
            }

            if (tail == capacity) {
                tail = 0;
            }

            array[tail] = item;
            size++;
            tail++;
            // Notify or notify all? There is only one consumer to profit from a just added item.
            lock.notify();
        }
    }

    @Override
    public final T take() throws InterruptedException {
        T item = null;

        synchronized (lock) {
            while (size == 0) {
                lock.wait();
            }

            if (head == capacity) {
                head = 0;
            }

            item = array[head];
            array[head] = null;
            head++;
            size--;

            lock.notifyAll();
        }

        return item;
    }

    @Override
    public boolean add(T t) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean offer(T t) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean offer(T t, long timeout, TimeUnit unit) throws InterruptedException {
        throw new UnsupportedOperationException("Not supported yet.");
    }


    @Override
    public T poll(long timeout, TimeUnit unit) throws InterruptedException {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int remainingCapacity() {
        return capacity - size;
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean contains(Object o) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int drainTo(Collection<? super T> c) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int drainTo(Collection<? super T> c, int maxElements) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public T remove() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public T poll() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public T element() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public T peek() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int size() {
        return size;
    }

    @Override
    public boolean isEmpty() {
        return size > 0;
    }

    @Override
    public Iterator<T> iterator() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object[] toArray() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public <T1> T1[] toArray(T1[] a) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean addAll(Collection<? extends T> c) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
