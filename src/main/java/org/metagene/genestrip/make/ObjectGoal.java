package org.metagene.genestrip.make;

public abstract class ObjectGoal<T,P> extends Goal<P> {
	private T object;
	
	@SafeVarargs
	public ObjectGoal(P project, String name, Goal<P>... dependencies) {
		super(project, name, dependencies);
	}
	
	public T get() {
		if (!isMade()) {
			make();
		}
		return object;
	}
	
	protected void set(T object) {
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Setting " + this + " to " + object);
		}
		this.object = object;
	}
	
	@Override
	public void cleanThis() {
		object = null;
	}

	@Override
	public boolean isMade() {
		return object != null;
	}
	
	@Override
	public String toString() {
		return "object goal: "+ getName();
	}
}
